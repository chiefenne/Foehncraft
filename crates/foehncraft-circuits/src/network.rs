//! The circuit model and its solution: the march of a tube in segments
//! and the Newton solve of the flows, pressures and enthalpies of the
//! whole graph.

use crate::case::{Boiling, Case, Chf, Friction, NodeKind, Subcooled, Void};
use crate::chf::{self, LookUpTable};
use crate::correlations as cor;
use crate::fluid::{self, Phase1, Point, SatProps};
use crate::{Error, G, codes, linalg};
use std::f64::consts::PI;

/// The largest specific enthalpy of a march (J/kg): steam near 1073 K at
/// low pressure; a tube with too little flow for its heat stops here.
const H_MAX: f64 = 4.1e6;

/// A node with its unknowns.
#[derive(Clone, Debug)]
pub struct NodeModel {
    /// Name.
    pub name: String,
    /// Kind.
    pub kind: NodeKind,
    /// Fixed pressure (Pa) of a pressure or drum node.
    pub p_fixed: Option<f64>,
    /// Enthalpy (J/kg) of the water a pressure node delivers, or of the
    /// inflow of a junction.
    pub h_given: Option<f64>,
    /// Mass flow (kg/s) entering at a junction.
    pub m_source: f64,
    /// Feedwater enthalpy (J/kg) of a drum.
    pub h_feed: f64,
    /// Index of the pressure unknown (junctions).
    pub ip: Option<usize>,
    /// Index of the enthalpy unknown (junctions and drums).
    pub ih: Option<usize>,
    /// Start pressure (Pa).
    pub p0: f64,
    /// Start enthalpy (J/kg).
    pub h0: f64,
}

/// A tube with its geometry and heat.
#[derive(Clone, Debug)]
pub struct TubeModel {
    /// Name.
    pub name: String,
    /// Node of positive flow entry.
    pub from: usize,
    /// Node of positive flow exit.
    pub to: usize,
    /// Length (m).
    pub length: f64,
    /// Inner diameter (m).
    pub d: f64,
    /// Flow area of one tube (m²).
    pub area: f64,
    /// Sine of the angle from `from` to `to` above the horizontal.
    pub sin: f64,
    /// Relative roughness.
    pub rel_rough: f64,
    /// Parallel tubes.
    pub count: f64,
    /// Inlet loss coefficient.
    pub loss: f64,
    /// Heat flux (W/m²) of each segment, ordered from `from` to `to`.
    pub q: Vec<f64>,
    /// Index of the mass flow unknown (kg/s of one tube).
    pub im: usize,
    /// Start flow of one tube (kg/s).
    pub m0: f64,
}

impl TubeModel {
    /// Segment length (m).
    pub fn dz(&self) -> f64 {
        self.length / self.q.len() as f64
    }

    /// Heat of one tube (W).
    pub fn heat(&self) -> f64 {
        self.q.iter().sum::<f64>() * PI * self.d * self.dz()
    }
}

/// The models of the case.
#[derive(Clone, Copy, Debug)]
pub struct Models {
    /// Two-phase friction.
    pub friction: Friction,
    /// Void.
    pub void: Void,
    /// Boiling heat transfer.
    pub boiling: Boiling,
    /// Subcooled void.
    pub subcooled: Subcooled,
    /// CHF source.
    pub chf: Chf,
}

/// The unknowns, the recorded marches and the Newton history.
pub type Solved = (Vec<f64>, Vec<March>, Vec<f64>);

/// A circuit ready to solve.
#[derive(Clone, Debug)]
pub struct Model {
    /// Nodes.
    pub nodes: Vec<NodeModel>,
    /// Tubes.
    pub tubes: Vec<TubeModel>,
    /// Number of unknowns.
    pub n: usize,
    /// Models.
    pub models: Models,
    /// The CHF table of `chf = "table"`.
    pub table: Option<LookUpTable>,
    /// Scale of the mass flows (kg/s).
    pub m_ref: f64,
    /// Newton tolerance.
    pub tolerance: f64,
    /// Newton iterations at most.
    pub max_iterations: usize,
    /// Smallest CHF ratio accepted without a warning.
    pub chf_margin: f64,
}

/// The local two-phase state at one point of a tube.
#[derive(Clone, Copy, Debug)]
pub struct Local {
    /// The state at (p, h).
    pub point: Point,
    /// Equilibrium quality (NaN above p_c).
    pub x_eq: f64,
    /// Flowing (true) quality.
    pub x: f64,
    /// Void fraction.
    pub alpha: f64,
    /// Mixture density for gravity (kg/m³).
    pub rho_mix: f64,
    /// Momentum density (kg/m³).
    pub rho_mom: f64,
    /// Friction pressure gradient (Pa/m), along the flow.
    pub dpdz_fric: f64,
    /// Saturated phases (below p_c).
    pub sat: Option<SatProps>,
    /// Liquid properties (the subcooled liquid in subcooled boiling).
    pub liquid: Phase1,
}

/// One segment of a solved tube, ordered along the flow.
#[derive(Clone, Debug)]
pub struct Segment {
    /// Distance from the `from` end to the middle of the segment (m).
    pub s: f64,
    /// Pressure at the middle (Pa).
    pub p: f64,
    /// Enthalpy at the middle (J/kg).
    pub h: f64,
    /// Bulk temperature (K).
    pub t: f64,
    /// Equilibrium quality.
    pub x_eq: f64,
    /// Flowing quality.
    pub x: f64,
    /// Void fraction.
    pub alpha: f64,
    /// Heat flux (W/m²).
    pub q: f64,
    /// Inner wall temperature (K).
    pub t_wall: f64,
    /// Heat transfer coefficient q/(T_w − T_b) (W/(m² K)), NaN unheated.
    pub htc: f64,
    /// Heat transfer regime.
    pub regime: &'static str,
    /// Critical heat flux (W/m²), NaN where not computed.
    pub q_chf: f64,
    /// CHF ratio q_CHF/q, NaN where not computed.
    pub chf_ratio: f64,
}

/// The result of the march of one tube.
#[derive(Clone, Debug, Default)]
pub struct March {
    /// Pressure at the downstream end (Pa).
    pub p_end: f64,
    /// Enthalpy at the downstream end (J/kg).
    pub h_end: f64,
    /// Pressure drop by friction (Pa), loss included.
    pub dp_fric: f64,
    /// Pressure drop by gravity (Pa).
    pub dp_grav: f64,
    /// Pressure drop by acceleration (Pa).
    pub dp_acc: f64,
    /// The segments (when recorded).
    pub segments: Vec<Segment>,
    /// The march hit the enthalpy limit.
    pub clipped: bool,
}

impl Model {
    /// Build the model of a case (reads the CHF table if one is named).
    pub fn new(case: &Case) -> Result<Self, Error> {
        let bad = |m: String| Error::new(codes::CASE, m);
        let enthalpy =
            |p: f64, t: Option<f64>, h: Option<f64>, what: &str| -> Result<Option<f64>, Error> {
                if let Some(h) = h {
                    return Ok(Some(h));
                }
                match t {
                    None => Ok(None),
                    Some(t) => foehncraft_kernels::steam::state_pt(p, t)
                        .map(|s| Some(s.h))
                        .ok_or_else(|| {
                            bad(format!(
                                "{what}: no IAPWS-IF97 state at p = {p} Pa, T = {t} K"
                            ))
                        }),
                }
            };
        let fixed: Vec<f64> = case
            .nodes
            .iter()
            .filter_map(|n| match n.kind {
                NodeKind::Junction => None,
                _ => n.pressure,
            })
            .collect();
        let p_mean = fixed.iter().sum::<f64>() / fixed.len() as f64;
        let mut nodes = Vec::new();
        let mut n = 0;
        // the start enthalpy: the mean of the given enthalpies
        let mut given = Vec::new();
        for c in &case.nodes {
            let p = c.pressure.unwrap_or(p_mean);
            let h = enthalpy(p, c.temperature, c.enthalpy, &format!("node `{}`", c.name))?;
            let h_feed = match (c.kind, c.feed_temperature) {
                (NodeKind::Drum, Some(t)) => {
                    enthalpy(p, Some(t), None, &format!("drum `{}`", c.name))?.unwrap_or(0.0)
                }
                _ => 0.0,
            };
            if let Some(h) = h {
                given.push(h);
            }
            if c.kind == NodeKind::Drum {
                let sat = foehncraft_kernels::steam::saturation_p(p).ok_or_else(|| {
                    bad(format!("drum `{}`: no saturation line at {p} Pa", c.name))
                })?;
                given.push(sat.liquid.h);
            }
            nodes.push(NodeModel {
                name: c.name.clone(),
                kind: c.kind,
                p_fixed: match c.kind {
                    NodeKind::Junction => None,
                    _ => c.pressure,
                },
                h_given: h,
                m_source: c.mass_flow.unwrap_or(0.0),
                h_feed,
                ip: None,
                ih: None,
                p0: p,
                h0: h.unwrap_or(f64::NAN),
            });
        }
        let h_mean = if given.is_empty() {
            1e6
        } else {
            given.iter().sum::<f64>() / given.len() as f64
        };
        let index = |name: &str| {
            nodes
                .iter()
                .position(|x: &NodeModel| x.name == name)
                .unwrap_or(0)
        };
        let mut tubes = Vec::new();
        for t in &case.tubes {
            let segs = t.segments.unwrap_or(case.circuit.segments);
            let dz = t.length / segs as f64;
            let area = PI * t.diameter * t.diameter / 4.0;
            let mut q = vec![0.0; segs];
            let heated = t.length - t.unheated_inlet;
            let shape = |s: f64| -> f64 {
                match &t.heat_shape {
                    None => 1.0,
                    Some(v) => {
                        let i = v.partition_point(|a| a[0] <= s).clamp(1, v.len() - 1) - 1;
                        let w = (s - v[i][0]) / (v[i + 1][0] - v[i][0]);
                        v[i][1] + w * (v[i + 1][1] - v[i][1])
                    }
                }
            };
            for (i, qi) in q.iter_mut().enumerate() {
                let z = (i as f64 + 0.5) * dz;
                if z >= t.unheated_inlet {
                    *qi = shape(z / t.length);
                }
            }
            let sum: f64 = q.iter().sum::<f64>() * dz;
            let q_mean = match (t.heat_flux, t.heat) {
                (Some(f), _) => f,
                (None, Some(w)) => w / (PI * t.diameter * heated),
                _ => 0.0,
            };
            if sum > 0.0 {
                let scale = q_mean * heated / sum;
                q.iter_mut().for_each(|v| *v *= scale);
            }
            tubes.push(TubeModel {
                name: t.name.clone(),
                from: index(&t.from),
                to: index(&t.to),
                length: t.length,
                d: t.diameter,
                area,
                sin: t.rise / t.length,
                rel_rough: t.roughness / t.diameter,
                count: t.count as f64,
                loss: t.loss,
                q,
                im: n,
                m0: t.start_flow.unwrap_or(1000.0 * area),
            });
            n += 1;
        }
        for nd in nodes.iter_mut() {
            match nd.kind {
                NodeKind::Junction => {
                    nd.ip = Some(n);
                    nd.ih = Some(n + 1);
                    n += 2;
                    if nd.h0.is_nan() {
                        nd.h0 = h_mean;
                    }
                }
                NodeKind::Drum => {
                    nd.ih = Some(n);
                    n += 1;
                    nd.h0 = foehncraft_kernels::steam::saturation_p(nd.p0)
                        .map_or(h_mean, |s| s.liquid.h);
                }
                NodeKind::Pressure => {}
            }
        }
        // every junction must be reached by a tube
        for (i, nd) in nodes.iter().enumerate() {
            if !tubes.iter().any(|t| t.from == i || t.to == i) {
                return Err(bad(format!("node `{}` has no tube", nd.name)));
            }
        }
        let table = match case.circuit.chf {
            Chf::Table => {
                let path = case.chf_table_path().unwrap_or_default();
                let text = std::fs::read_to_string(&path)
                    .map_err(|e| bad(format!("CHF table {}: {e}", path.display())))?;
                Some(LookUpTable::from_csv(&text)?)
            }
            _ => None,
        };
        let src: f64 = nodes.iter().map(|x| x.m_source.abs()).sum();
        let flows: f64 =
            tubes.iter().map(|t| t.m0.abs() * t.count).sum::<f64>() / tubes.len() as f64;
        Ok(Self {
            nodes,
            tubes,
            n,
            models: Models {
                friction: case.circuit.friction,
                void: case.circuit.void,
                boiling: case.circuit.boiling,
                subcooled: case.circuit.subcooled,
                chf: case.circuit.chf,
            },
            table,
            m_ref: src.max(flows).max(1e-6),
            tolerance: case.circuit.tolerance,
            max_iterations: case.circuit.max_iterations,
            chf_margin: case.circuit.chf_margin,
        })
    }

    /// The start vector of the unknowns.
    pub fn start(&self) -> Vec<f64> {
        let mut x = vec![0.0; self.n];
        for t in &self.tubes {
            x[t.im] = t.m0;
        }
        for nd in &self.nodes {
            if let Some(i) = nd.ip {
                x[i] = nd.p0;
            }
            if let Some(i) = nd.ih {
                x[i] = nd.h0;
            }
        }
        x
    }

    /// Pressure of a node for the unknowns x.
    pub fn p(&self, x: &[f64], node: usize) -> f64 {
        let nd = &self.nodes[node];
        nd.p_fixed.unwrap_or_else(|| x[nd.ip.unwrap_or(0)])
    }

    /// Enthalpy of the water leaving a node into a tube.
    pub fn h(&self, x: &[f64], node: usize) -> f64 {
        let nd = &self.nodes[node];
        match nd.ih {
            Some(i) => x[i],
            None => nd.h_given.unwrap_or(f64::NAN),
        }
    }

    /// The local state at (p, h) with mass flux g, heat flux q and the
    /// sine `sin_up` of the flow direction.
    pub fn local(
        &self,
        t: &TubeModel,
        p: f64,
        h: f64,
        g: f64,
        q: f64,
        sin_up: f64,
    ) -> Result<Local, Error> {
        let point = Point::new(p, h)?;
        let rough = t.rel_rough;
        let single = |s: &foehncraft_kernels::steam::State| -> Local {
            let pr = Phase1::of(s);
            let re = g * t.d / pr.mu;
            let f = cor::darcy_churchill(re, rough);
            Local {
                point,
                x_eq: point.x_eq(),
                x: if point.x_eq() >= 1.0 { 1.0 } else { 0.0 },
                alpha: if point.x_eq() >= 1.0 { 1.0 } else { 0.0 },
                rho_mix: s.rho,
                rho_mom: s.rho,
                dpdz_fric: f * g * g / (2.0 * t.d * s.rho),
                sat: point.sat.as_ref().map(SatProps::of),
                liquid: pr,
            }
        };
        let Some(sat) = point.sat else {
            return Ok(single(&point.state));
        };
        let x_eq = point.x_eq();
        if x_eq >= 1.0 {
            return Ok(single(&point.state));
        }
        let sp = SatProps::of(&sat);
        let x = match self.models.subcooled {
            Subcooled::SahaZuber if q > 0.0 && x_eq < 1.0 => {
                let xd = cor::saha_zuber_x_d(q, g, t.d, sp.l.k, sp.l.cp, sp.h_fg);
                cor::levy_quality(x_eq, xd).min(1.0)
            }
            _ => x_eq.clamp(0.0, 1.0),
        };
        if x <= 0.0 {
            let mut l = single(&point.state);
            l.x_eq = x_eq;
            l.sat = Some(sp);
            return Ok(l);
        }
        // the liquid of subcooled boiling is subcooled: h_l = (h − x h″)/(1 − x)
        let liquid = if x_eq < 0.0 {
            let hl = (h - x * sat.vapour.h) / (1.0 - x);
            match foehncraft_kernels::steam::state_ph(p, hl.min(sat.liquid.h)) {
                Some(s) if s.cp.is_finite() => Phase1::of(&s),
                _ => sp.l,
            }
        } else {
            sp.l
        };
        let (rl, rg) = (liquid.rho, sp.g.rho);
        let alpha = cor::void_fraction(
            x,
            g,
            rl,
            rg,
            sp.sigma,
            sin_up,
            self.models.void == Void::DriftFlux,
        );
        let rho_mix = alpha * rg + (1.0 - alpha) * rl;
        let rho_mom = cor::momentum_density(x, alpha, rl, rg);
        let dpdz_fric = match self.models.friction {
            Friction::Friedel => {
                let f_lo = cor::darcy_churchill(g * t.d / liquid.mu, rough);
                let f_go = cor::darcy_churchill(g * t.d / sp.g.mu, rough);
                let phi2 = cor::friedel_multiplier(
                    x, g, t.d, rl, rg, liquid.mu, sp.g.mu, sp.sigma, f_lo, f_go,
                );
                phi2 * f_lo * g * g / (2.0 * t.d * rl)
            }
            Friction::Homogeneous => {
                let mu = 1.0 / (x / sp.g.mu + (1.0 - x) / liquid.mu);
                let rho_h = 1.0 / (x / rg + (1.0 - x) / rl);
                cor::darcy_churchill(g * t.d / mu, rough) * g * g / (2.0 * t.d * rho_h)
            }
        };
        Ok(Local {
            point,
            x_eq,
            x,
            alpha,
            rho_mix,
            rho_mom,
            dpdz_fric,
            sat: Some(sp),
            liquid,
        })
    }

    /// March one tube with the flow m of one tube (kg/s; negative from
    /// `to` to `from`) from the upstream pressure and enthalpy.
    pub fn march(
        &self,
        t: &TubeModel,
        m: f64,
        p_up: f64,
        h_up: f64,
        record: bool,
    ) -> Result<March, Error> {
        let forward = m >= 0.0;
        let ma = m.abs().max(1e-12);
        let g = ma / t.area;
        let sin_up = if forward { t.sin } else { -t.sin };
        let nseg = t.q.len();
        let dz = t.dz();
        let mut out = March::default();
        let (mut p, mut h) = (p_up, h_up);
        let mut here = self.local(t, p, h, g, t.q[if forward { 0 } else { nseg - 1 }], sin_up)?;
        if t.loss > 0.0 {
            let dp = t.loss * g * g / (2.0 * here.rho_mix);
            p -= dp;
            out.dp_fric += dp;
            here = self.local(t, p, h, g, t.q[if forward { 0 } else { nseg - 1 }], sin_up)?;
        }
        let mut dp_last = 0.0;
        for k in 0..nseg {
            let i = if forward { k } else { nseg - 1 - k };
            let q = t.q[i];
            let mut h_out = h + q * PI * t.d * dz / ma;
            if h_out > H_MAX {
                h_out = H_MAX;
                out.clipped = true;
            }
            let mut p_out = p - dp_last;
            let mut there = here;
            let (mut df, mut dg, mut da) = (0.0, 0.0, 0.0);
            for _ in 0..3 {
                there = self.local(t, p_out, h_out, g, q, sin_up)?;
                df = 0.5 * (here.dpdz_fric + there.dpdz_fric) * dz;
                dg = 0.5 * (here.rho_mix + there.rho_mix) * G * sin_up * dz;
                da = g * g * (1.0 / there.rho_mom - 1.0 / here.rho_mom);
                p_out = p - df - dg - da;
                if p_out <= 1e3 {
                    return Err(Error::new(
                        codes::RANGE,
                        format!(
                            "tube `{}`: the pressure falls below 1 kPa (too much flow for the pressure level)",
                            t.name
                        ),
                    ));
                }
            }
            out.dp_fric += df;
            out.dp_grav += dg;
            out.dp_acc += da;
            dp_last = p - p_out;
            if record {
                let s_mid = (i as f64 + 0.5) * dz;
                let mid = self.local(t, 0.5 * (p + p_out), 0.5 * (h + h_out), g, q, sin_up)?;
                out.segments.push(self.segment(t, &mid, g, q, s_mid));
            }
            p = p_out;
            h = h_out;
            here = there;
        }
        out.p_end = p;
        out.h_end = h;
        Ok(out)
    }

    /// The wall of segment `i` of tube `ti` at the pressure p and enthalpy
    /// h of its middle with the flow m of one tube (kg/s), for the heat
    /// flux q (W/m²) in place of the segment's own: the heat transfer
    /// coefficient of a wall coupled to a 3D mesh (SF7c) where the heat
    /// of the last solve was zero or small.
    pub fn wall_at(
        &self,
        ti: usize,
        i: usize,
        p: f64,
        h: f64,
        m: f64,
        q: f64,
    ) -> Result<Segment, Error> {
        let t = &self.tubes[ti];
        let g = m.abs().max(1e-12) / t.area;
        let sin_up = if m >= 0.0 { t.sin } else { -t.sin };
        let l = self.local(t, p, h, g, q, sin_up)?;
        Ok(self.segment(t, &l, g, q, (i as f64 + 0.5) * t.dz()))
    }

    /// The wall temperature, regime and CHF of a segment.
    fn segment(&self, t: &TubeModel, l: &Local, g: f64, q: f64, s: f64) -> Segment {
        let tb = l.point.t();
        let mut seg = Segment {
            s,
            p: l.point.p,
            h: l.point.h,
            t: tb,
            x_eq: l.x_eq,
            x: l.x,
            alpha: l.alpha,
            q,
            t_wall: tb,
            htc: f64::NAN,
            regime: "unheated",
            q_chf: f64::NAN,
            chf_ratio: f64::NAN,
        };
        if q <= 0.0 {
            return seg;
        }
        let gnielinski = |ph: &Phase1, re: f64| cor::nusselt_tube(re, ph.pr()) * ph.k / t.d;
        let sp = match &l.sat {
            None => {
                let ph = Phase1::of(&l.point.state);
                let hc = gnielinski(&ph, g * t.d / ph.mu);
                seg.t_wall = tb + q / hc;
                seg.regime = "supercritical";
                seg.htc = hc;
                return seg;
            }
            Some(sp) => *sp,
        };
        if l.x_eq >= 1.0 {
            let ph = Phase1::of(&l.point.state);
            let hc = gnielinski(&ph, g * t.d / ph.mu);
            seg.t_wall = tb + q / hc;
            seg.regime = "superheated vapour";
            seg.htc = hc;
            return seg;
        }
        // the critical heat flux at the local conditions
        let q_chf = match self.models.chf {
            Chf::Bowring => chf::bowring(l.point.p, g, t.d, sp.h_fg, l.x_eq),
            Chf::Table => self
                .table
                .as_ref()
                .map_or(f64::NAN, |tb| tb.value(l.point.p, g, l.x_eq, t.d)),
            Chf::None => f64::NAN,
        };
        seg.q_chf = q_chf;
        seg.chf_ratio = q_chf / q;
        let input = |x: f64, t_bulk: f64| cor::BoilingInput {
            q,
            g_flux: g,
            d: t.d,
            x,
            t_bulk,
            p: l.point.p,
            rho_l: sp.l.rho,
            rho_g: sp.g.rho,
            mu_l: sp.l.mu,
            mu_g: sp.g.mu,
            k_l: sp.l.k,
            cp_l: sp.l.cp,
            sigma: sp.sigma,
            h_fg: sp.h_fg,
            t_sat: sp.t,
        };
        let boil = |b: &cor::BoilingInput| match self.models.boiling {
            Boiling::Chen => cor::chen_wall_temperature(b, fluid::p_sat),
            Boiling::GungorWinterton => {
                cor::gungor_winterton_wall_temperature(b, fluid::P_CRIT, t.sin.abs() < 0.5)
            }
        };
        if seg.chf_ratio < 1.0 {
            // beyond the CHF: vapour convection of Dougall & Rohsenow (1963)
            let x = l.x_eq.clamp(0.0, 1.0);
            let re = g * t.d / sp.g.mu * (x + sp.g.rho / sp.l.rho * (1.0 - x));
            let hc = cor::dittus_boelter(re, sp.g.pr(), sp.g.k, t.d);
            seg.t_wall = sp.t + q / hc;
            seg.htc = q / (seg.t_wall - tb);
            seg.regime = "post-CHF";
            return seg;
        }
        if l.x_eq < 0.0 {
            let ph = Phase1::of(&l.point.state);
            let hc = gnielinski(&ph, g * t.d / ph.mu);
            let tw = tb + q / hc;
            if tw <= sp.t {
                seg.t_wall = tw;
                seg.htc = hc;
                seg.regime = "liquid";
            } else {
                seg.t_wall = boil(&input(l.x_eq, tb)).max(sp.t).min(tw);
                seg.htc = q / (seg.t_wall - tb);
                seg.regime = "subcooled boiling";
            }
            return seg;
        }
        seg.t_wall = boil(&input(l.x_eq, sp.t));
        seg.htc = q / (seg.t_wall - tb).max(1e-9);
        seg.regime = "saturated boiling";
        seg
    }

    /// Upstream state of a tube for the unknowns x: (pressure, enthalpy,
    /// downstream node).
    fn upstream(&self, x: &[f64], t: &TubeModel) -> (f64, f64, usize) {
        let m = x[t.im];
        if m >= 0.0 {
            (self.p(x, t.from), self.h(x, t.from), t.to)
        } else {
            (self.p(x, t.to), self.h(x, t.to), t.from)
        }
    }

    /// The scaled residuals at x, and the marches.
    pub fn residual(&self, x: &[f64], record: bool) -> Result<(Vec<f64>, Vec<March>), Error> {
        const P_SCALE: f64 = 1e5;
        const H_SCALE: f64 = 1e5;
        let mut r = vec![0.0; self.n];
        let mut marches = Vec::with_capacity(self.tubes.len());
        for t in &self.tubes {
            let (p_up, h_up, down) = self.upstream(x, t);
            let mr = self.march(t, x[t.im], p_up, h_up, record)?;
            r[t.im] = (mr.p_end - self.p(x, down)) / P_SCALE;
            marches.push(mr);
        }
        // node balances
        let mref = self.m_ref;
        for (k, nd) in self.nodes.iter().enumerate() {
            // inflows: (mass flow, arriving enthalpy); outflow total
            let mut m_net = nd.m_source;
            let mut inflow = Vec::new();
            if nd.m_source > 0.0 {
                inflow.push((nd.m_source, nd.h_given.unwrap_or(f64::NAN)));
            }
            for (t, mr) in self.tubes.iter().zip(&marches) {
                let m = x[t.im] * t.count;
                if t.to == k {
                    m_net += m;
                    if m > 0.0 {
                        inflow.push((m, mr.h_end));
                    }
                }
                if t.from == k {
                    m_net -= m;
                    if m < 0.0 {
                        inflow.push((-m, mr.h_end));
                    }
                }
            }
            match nd.kind {
                NodeKind::Junction => {
                    let (ip, ih) = (nd.ip.unwrap_or(0), nd.ih.unwrap_or(0));
                    r[ip] = m_net / mref;
                    let hn = x[ih];
                    let e: f64 = inflow.iter().map(|(m, h)| m * (h - hn)).sum();
                    r[ih] = e / (mref * H_SCALE);
                }
                NodeKind::Drum => {
                    let ih = nd.ih.unwrap_or(0);
                    let hn = x[ih];
                    let sat = foehncraft_kernels::steam::saturation_p(nd.p_fixed.unwrap_or(0.0));
                    let (hf, hg) = sat.map_or((hn, hn), |s| (s.liquid.h, s.vapour.h));
                    let (mut m_liq, mut e_liq, mut m_steam) = (0.0, 0.0, 0.0);
                    for (m, h) in &inflow {
                        let xq = ((h - hf) / (hg - hf)).clamp(0.0, 1.0);
                        m_steam += m * xq;
                        m_liq += m * (1.0 - xq);
                        e_liq += m * (1.0 - xq) * if xq > 0.0 { hf } else { *h };
                    }
                    let e = e_liq + m_steam * nd.h_feed - hn * (m_liq + m_steam);
                    r[ih] = e / (mref * H_SCALE);
                }
                NodeKind::Pressure => {}
            }
        }
        Ok((r, marches))
    }

    /// Solve the circuit by a damped Newton method with a finite
    /// difference Jacobian; returns the unknowns, the marches (recorded)
    /// and the Newton history (largest scaled residual per iteration).
    pub fn solve(&self) -> Result<Solved, Error> {
        self.solve_from(self.start())
    }

    /// [`Model::solve`] from the unknowns `x` (a warm start: the solution
    /// of a circuit whose heat changed a little, SF7c).
    pub fn solve_from(&self, mut x: Vec<f64>) -> Result<Solved, Error> {
        let norm = |r: &[f64]| r.iter().fold(0.0_f64, |m, v| m.max(v.abs()));
        let (mut r, _) = self.residual(&x, false)?;
        let mut history = vec![norm(&r)];
        let n = self.n;
        let typical: Vec<f64> = x.iter().map(|v| v.abs().max(1e-6)).collect();
        for _ in 0..self.max_iterations {
            if norm(&r) <= self.tolerance {
                let (_, marches) = self.residual(&x, true)?;
                return Ok((x, marches, history));
            }
            let mut jac = vec![0.0; n * n];
            for j in 0..n {
                let step = 1e-7 * x[j].abs().max(typical[j]);
                let mut xp = x.clone();
                xp[j] += step;
                let (rp, _) = self.residual(&xp, false)?;
                for i in 0..n {
                    jac[i * n + j] = (rp[i] - r[i]) / step;
                }
            }
            let mut dx: Vec<f64> = r.iter().map(|v| -v).collect();
            linalg::solve(&mut jac, &mut dx, n).ok_or_else(|| {
                Error::new(
                    codes::CONVERGENCE,
                    "the Jacobian of the circuit is singular: a part of the circuit is not connected to a pressure or drum node, or a tube carries no flow",
                )
            })?;
            // limit the step of the enthalpies and the pressures
            let mut lam: f64 = 1.0;
            for nd in &self.nodes {
                if let Some(i) = nd.ih {
                    lam = lam.min(3e5 / dx[i].abs().max(1e-300));
                }
                if let Some(i) = nd.ip {
                    lam = lam.min(0.3 * x[i] / dx[i].abs().max(1e-300));
                }
            }
            let r0 = norm(&r);
            let mut accepted = false;
            for _ in 0..30 {
                let xn: Vec<f64> = x.iter().zip(&dx).map(|(a, d)| a + lam * d).collect();
                if let Ok((rn, _)) = self.residual(&xn, false)
                    && (norm(&rn) < r0 * (1.0 - 1e-4 * lam) || norm(&rn) <= self.tolerance)
                {
                    x = xn;
                    r = rn;
                    accepted = true;
                    break;
                }
                lam *= 0.5;
            }
            if !accepted {
                return Err(Error::new(
                    codes::CONVERGENCE,
                    format!(
                        "the Newton solve of the circuit stalled at a residual of {r0:.3e} (give start flows with `start_flow`, or check that the heat can be carried by the flow)"
                    ),
                ));
            }
            history.push(norm(&r));
        }
        Err(Error::new(
            codes::CONVERGENCE,
            format!(
                "the circuit did not converge in {} Newton iterations (residual {:.3e}, tolerance {:.1e})",
                self.max_iterations,
                norm(&r),
                self.tolerance
            ),
        ))
    }
}
