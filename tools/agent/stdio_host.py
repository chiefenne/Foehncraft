"""A minimal MCP host for scripts and for a coding assistant that has no MCP
registration of its own: it starts `foehn mcp`, keeps it open, and lets any
number of short commands send tool calls to it through a Unix socket, even
while a long `flow.run` is going (the control port, L19). Every call and
answer goes into a transcript (JSON lines).

    python3 -m tools.agent.stdio_host serve SOCKET TRANSCRIPT -- foehn --quiet mcp ...
    python3 -m tools.agent.stdio_host call SOCKET TOOL '{"json": "arguments"}'
    python3 -m tools.agent.stdio_host start SOCKET TOOL '{"json": "arguments"}'   # prints a ticket
    python3 -m tools.agent.stdio_host wait SOCKET TICKET
    python3 -m tools.agent.stdio_host stop SOCKET

The host announces no capabilities, so a call that needs the person's
approval (`ask` in the action table) is refused (`MCP005`); the person
allows what a study may do in `foehn-mcp.toml` beforehand. `call` prints the
first line of the answer and the result as JSON; `start` returns at once
(for a call that takes minutes) and `wait` collects its answer.
"""
import json
import os
import socket
import subprocess
import sys
import threading
import time


def serve(sock_path, transcript, argv):
    p = subprocess.Popen(argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1)
    lock = threading.Lock()
    waiting = {}
    state = {"n": 0}
    log = open(transcript, "a", buffering=1)

    def write(m):
        with lock:
            p.stdin.write(json.dumps(m) + "\n")
            p.stdin.flush()

    def reader():
        for line in p.stdout:
            m = json.loads(line)
            if "id" in m and "method" not in m and m["id"] in waiting:
                waiting[m["id"]]["answer"] = m
                waiting[m["id"]]["event"].set()

    threading.Thread(target=reader, daemon=True).start()

    def post(method, params):
        with lock:
            state["n"] += 1
            i = state["n"]
        w = {"event": threading.Event(), "answer": None}
        waiting[i] = w
        write({"jsonrpc": "2.0", "id": i, "method": method, "params": params})
        return i, w

    post_i = post("initialize", {"protocolVersion": "2025-06-18", "capabilities": {},
                                 "clientInfo": {"name": "stdio-host", "version": "1"}})
    post_i[1]["event"].wait()
    write({"jsonrpc": "2.0", "method": "notifications/initialized"})
    tickets = {}

    def finish(ticket, tool, args, i, w, t0):
        w["event"].wait()
        r = w["answer"]["result"]
        env = r.get("structuredContent")
        out = {"ticket": ticket, "tool": tool, "arguments": args,
               "started": t0, "seconds": round(time.time() - t0, 2),
               "headline": r["content"][0]["text"].split("\n")[0], "envelope": env}
        log.write(json.dumps(out) + "\n")
        return out

    def handle(conn):
        req = json.loads(conn.makefile().readline())
        op = req["op"]
        if op == "call" or op == "start":
            with lock:
                ticket = f"t{len(tickets) + 1}"
                tickets[ticket] = None
            t0 = time.time()
            i, w = post("tools/call", {"name": req["tool"], "arguments": req["arguments"]})
            if op == "start":
                def bg():
                    tickets[ticket] = finish(ticket, req["tool"], req["arguments"], i, w, t0)
                threading.Thread(target=bg, daemon=True).start()
                conn.sendall((json.dumps({"ticket": ticket}) + "\n").encode())
            else:
                out = finish(ticket, req["tool"], req["arguments"], i, w, t0)
                conn.sendall((json.dumps(out) + "\n").encode())
        elif op == "tools-list":
            i, w = post("tools/list", {})
            w["event"].wait()
            conn.sendall((json.dumps(w["answer"]["result"]) + "\n").encode())
        elif op == "wait":
            while tickets.get(req["ticket"]) is None:
                time.sleep(0.2)
            conn.sendall((json.dumps(tickets[req["ticket"]]) + "\n").encode())
        elif op == "stop":
            conn.sendall(b'{"ok": true}\n')
            p.stdin.close()
            p.wait()
            os._exit(0)
        conn.close()

    if os.path.exists(sock_path):
        os.remove(sock_path)
    s = socket.socket(socket.AF_UNIX)
    s.bind(sock_path)
    s.listen(8)
    while True:
        c, _ = s.accept()
        threading.Thread(target=handle, args=(c,), daemon=True).start()


def client(sock_path, req):
    s = socket.socket(socket.AF_UNIX)
    s.connect(sock_path)
    s.sendall((json.dumps(req) + "\n").encode())
    data = b""
    while not data.endswith(b"\n"):
        chunk = s.recv(1 << 16)
        if not chunk:
            break
        data += chunk
    return json.loads(data)


def show(out, full=False):
    print(out["headline"])
    env = out["envelope"] or {}
    if full:
        print(json.dumps(env.get("result"), indent=1))
    else:
        res = json.dumps(env.get("result"))
        print(res if len(res) < 2400 else res[:2400] + " ...")
    for d in env.get("diagnostics", []):
        print(f"  {d['code']}: {d['message']}")


def main():
    a = sys.argv[1:]
    if not a:
        sys.exit(__doc__)
    if a[0] == "serve":
        i = a.index("--")
        serve(a[1], a[2], a[i + 1:])
    elif a[0] in ("call", "start"):
        out = client(a[1], {"op": a[0], "tool": a[2], "arguments": json.loads(a[3]) if len(a) > 3 else {}})
        if a[0] == "start":
            print(out["ticket"])
        else:
            show(out, "--full" in a)
    elif a[0] == "wait":
        show(client(a[1], {"op": "wait", "ticket": a[2]}))
    elif a[0] == "tools":
        print(" ".join(t["name"] for t in client(a[1], {"op": "tools-list"})["tools"]))
    elif a[0] == "stop":
        client(a[1], {"op": "stop"})


if __name__ == "__main__":
    main()
