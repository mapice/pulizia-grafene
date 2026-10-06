"""Bound own scientific child process; never touches other user's processes.

macOS RSS is checked throughout initialization and calculation. Child gets
its own process group; only that group is terminated on cap or interruption.
"""
from pathlib import Path
import os,subprocess,time,signal,argparse,json

def main():
    p=argparse.ArgumentParser();p.add_argument('--mb',type=int,default=3000);p.add_argument('--minutes',type=float,default=30);p.add_argument('command',nargs=argparse.REMAINDER)
    args=p.parse_args();command=args.command
    if command and command[0]=='--':command=command[1:]
    if not command:raise SystemExit('Missing child command')
    child=subprocess.Popen(command,start_new_session=True)
    start=time.monotonic();peak=0;reason=None
    try:
        while child.poll() is None:
            r=subprocess.run(['ps','-o','rss=','-p',str(child.pid)],capture_output=True,text=True)
            if r.stdout.strip():peak=max(peak,float(r.stdout.strip())/1024)
            if peak>args.mb:reason='RSS cap'
            if time.monotonic()-start>args.minutes*60:reason='Wall-time cap'
            if reason:
                os.killpg(child.pid,signal.SIGTERM)
                try:child.wait(timeout=10)
                except subprocess.TimeoutExpired:os.killpg(child.pid,signal.SIGKILL);child.wait()
                break
            time.sleep(2)
    except BaseException:
        if child.poll() is None:os.killpg(child.pid,signal.SIGTERM)
        child.wait();raise
    print(json.dumps(dict(guard=dict(command=command,pid=child.pid,exitcode=child.returncode,peak_RSS_MB=peak,elapsed_seconds=time.monotonic()-start,termination_reason=reason))),flush=True)
    code=(128+abs(child.returncode)) if child.returncode and child.returncode<0 else child.returncode
    raise SystemExit(code if code else 1 if reason else 0)

if __name__=='__main__':main()
