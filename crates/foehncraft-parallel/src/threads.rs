//! Ranks as threads of one process: a [`Communicator`] for tests and tools.
//!
//! [`ThreadComm::run`] starts one thread per rank and hands each its
//! communicator. Messages travel through one unbounded channel per ordered
//! pair of ranks, so they arrive in the order they were sent and a send
//! never blocks; collectives are built from those messages (gather on rank
//! 0, then broadcast), so every rank receives the same bits. A multi-rank
//! algorithm can therefore be tested bit for bit against one rank with
//! `cargo test`, without an MPI installation; CI runs the same algorithms
//! under MPI as well.

use crate::comm::Communicator;
use std::sync::mpsc::{Receiver, Sender, channel};
use std::sync::{Arc, Barrier, Mutex};

enum Msg {
    F64(Vec<f64>),
    I64(Vec<i64>),
}

/// The communicator of one rank of [`ThreadComm::run`].
pub struct ThreadComm {
    rank: usize,
    size: usize,
    /// `to[d]`: channel to rank `d`.
    to: Vec<Sender<Msg>>,
    /// `from[s]`: channel from rank `s` (behind a lock so the communicator
    /// can be shared, `Send + Sync`, like the MPI one).
    from: Vec<Mutex<Receiver<Msg>>>,
    barrier: Arc<Barrier>,
}

impl ThreadComm {
    /// Runs `f` on `size` ranks, one thread each, and returns the results
    /// in rank order. A panic on any rank is passed on.
    pub fn run<R: Send>(size: usize, f: impl Fn(Arc<ThreadComm>) -> R + Sync) -> Vec<R> {
        let f = &f;
        let comms = Self::comms(size);
        std::thread::scope(|s| {
            let handles: Vec<_> = comms
                .into_iter()
                .map(|c| {
                    std::thread::Builder::new()
                        .name(format!("rank {}", c.rank))
                        .spawn_scoped(s, move || f(Arc::new(c)))
                        .expect("spawn rank thread")
                })
                .collect();
            handles
                .into_iter()
                .map(|h| match h.join() {
                    Ok(r) => r,
                    Err(e) => std::panic::resume_unwind(e),
                })
                .collect()
        })
    }

    /// Runs `lead` as rank 0 on the calling thread (it need not be `Send`:
    /// it may hold the log and the run record of the caller) and `follow`
    /// on ranks 1 to `size - 1`, one thread each (MPI step W7). Returns the
    /// result of rank 0. When rank 0 ends early (an error), the other
    /// ranks end at their next message to it, and the other way round;
    /// the panics of ranks 1 and up are dropped, a panic of rank 0 is
    /// passed on.
    pub fn run_lead<R>(
        size: usize,
        lead: impl FnOnce(Arc<ThreadComm>) -> R,
        follow: impl Fn(Arc<ThreadComm>) + Sync,
    ) -> R {
        let follow = &follow;
        let mut comms = Self::comms(size).into_iter();
        let first = comms.next().expect("at least one rank");
        std::thread::scope(|s| {
            let handles: Vec<_> = comms
                .map(|c| {
                    std::thread::Builder::new()
                        .name(format!("rank {}", c.rank))
                        .spawn_scoped(s, move || follow(Arc::new(c)))
                        .expect("spawn rank thread")
                })
                .collect();
            let r =
                std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| lead(Arc::new(first))));
            // a rank that fails ends at its next message to rank 0, and
            // rank 0 at its next message to it: the run ends with the
            // result of rank 0, the panics of the others are dropped
            for h in handles {
                let _ = h.join();
            }
            match r {
                Ok(r) => r,
                Err(e) => std::panic::resume_unwind(e),
            }
        })
    }

    /// The communicators of `size` ranks, connected.
    fn comms(size: usize) -> Vec<ThreadComm> {
        assert!(size > 0, "at least one rank");
        // senders[s][d], receivers[d][s]
        let mut senders: Vec<Vec<Sender<Msg>>> = (0..size).map(|_| Vec::new()).collect();
        let mut receivers: Vec<Vec<Option<Receiver<Msg>>>> = (0..size)
            .map(|_| (0..size).map(|_| None).collect())
            .collect();
        for (s, row) in senders.iter_mut().enumerate() {
            for rx in receivers.iter_mut() {
                let (tx, r) = channel();
                row.push(tx);
                rx[s] = Some(r);
            }
        }
        let barrier = Arc::new(Barrier::new(size));
        let comms: Vec<ThreadComm> = senders
            .into_iter()
            .zip(receivers)
            .enumerate()
            .map(|(rank, (to, from))| ThreadComm {
                rank,
                size,
                to,
                from: from
                    .into_iter()
                    .map(|r| Mutex::new(r.expect("channel")))
                    .collect(),
                barrier: barrier.clone(),
            })
            .collect();
        comms
    }

    fn send(&self, dest: usize, m: Msg) {
        self.to[dest].send(m).expect("rank thread ended early");
    }

    fn receive(&self, source: usize) -> Msg {
        self.from[source]
            .lock()
            .unwrap_or_else(std::sync::PoisonError::into_inner)
            .recv()
            .expect("rank thread ended early")
    }

    fn recv_f64(&self, source: usize) -> Vec<f64> {
        match self.receive(source) {
            Msg::F64(v) => v,
            Msg::I64(_) => panic!("rank {}: expected f64 values from {source}", self.rank),
        }
    }

    fn recv_i64(&self, source: usize) -> Vec<i64> {
        match self.receive(source) {
            Msg::I64(v) => v,
            Msg::F64(_) => panic!("rank {}: expected i64 values from {source}", self.rank),
        }
    }
}

impl Communicator for ThreadComm {
    fn rank(&self) -> usize {
        self.rank
    }
    fn size(&self) -> usize {
        self.size
    }
    fn barrier(&self) {
        self.barrier.wait();
    }
    fn all_gather_f64(&self, local: f64) -> Vec<f64> {
        if self.rank == 0 {
            let mut all = vec![local];
            for s in 1..self.size {
                all.extend(self.recv_f64(s));
            }
            for d in 1..self.size {
                self.send(d, Msg::F64(all.clone()));
            }
            all
        } else {
            self.send(0, Msg::F64(vec![local]));
            self.recv_f64(0)
        }
    }
    fn send_recv_f64(&self, send: &[f64], dest: usize, recv: &mut [f64], source: usize) {
        self.send(dest, Msg::F64(send.to_vec()));
        let v = self.recv_f64(source);
        assert_eq!(v.len(), recv.len(), "message length from rank {source}");
        recv.copy_from_slice(&v);
    }
    fn all_sum_i64(&self, buf: &mut [i64]) {
        if self.rank == 0 {
            for s in 1..self.size {
                let v = self.recv_i64(s);
                for (b, x) in buf.iter_mut().zip(v) {
                    *b += x;
                }
            }
            for d in 1..self.size {
                self.send(d, Msg::I64(buf.to_vec()));
            }
        } else {
            self.send(0, Msg::I64(buf.to_vec()));
            buf.copy_from_slice(&self.recv_i64(0));
        }
    }
    fn name(&self) -> &'static str {
        "threads"
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use foehncraft_core::exact::ExactSum;

    #[test]
    fn lead_on_the_calling_thread() {
        // the lead may hold what is not Send (a log on the caller's stack)
        let log = std::rc::Rc::new(std::cell::RefCell::new(Vec::new()));
        let l = log.clone();
        let sum = ThreadComm::run_lead(
            3,
            move |c| {
                let s = c.max(0.0) + c.all_gather_f64(1.0).iter().sum::<f64>();
                l.borrow_mut().push(s);
                s
            },
            |c| {
                c.max(c.rank() as f64);
                c.all_gather_f64(1.0);
            },
        );
        assert_eq!(sum, 5.0);
        assert_eq!(log.borrow().as_slice(), &[5.0]);
        // rank 0 ends early: the others end too
        let r = ThreadComm::run_lead(
            3,
            |_| 7,
            |c| {
                c.max(1.0);
            },
        );
        assert_eq!(r, 7);
    }

    #[test]
    fn ranks_exchange_and_reduce() {
        let out = ThreadComm::run(4, |c| {
            let r = c.rank();
            // ring
            let mut got = [0.0; 2];
            let next = (r + 1) % c.size();
            let prev = (r + c.size() - 1) % c.size();
            c.send_recv_f64(&[r as f64, 10.0 * r as f64], next, &mut got, prev);
            assert_eq!(got, [prev as f64, 10.0 * prev as f64]);
            let g = c.all_gather_f64(r as f64);
            assert_eq!(g, [0.0, 1.0, 2.0, 3.0]);
            let mut acc = ExactSum::new();
            acc.add(if r == 0 {
                1e300
            } else if r == 3 {
                -1e300
            } else {
                0.5
            });
            c.barrier();
            let local: Vec<f64> = (0..r).map(|i| (10 * r + i) as f64).collect();
            let all = c.all_gather_vec_f64(&local);
            for (q, v) in all.iter().enumerate() {
                assert_eq!(v, &(0..q).map(|i| (10 * q + i) as f64).collect::<Vec<_>>());
            }
            (c.sum_exact(&acc), c.max(r as f64), c.min(r as f64))
        });
        assert!(out.iter().all(|&x| x == (1.0, 3.0, 0.0)));
    }
}
