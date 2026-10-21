"""GPU-queue gate: exit 0 (and print OK) only when M35_DONE and the user-release flag GO_STAGE1 exist AND no M32/M33/M35 process is left on this machine (a process counts when its command
line or its working directory mentions M32-C1-PLANTCAD2 / M33-NTV2-FOLLOWUP; transient probes/waiters and this gate itself are ignored) AND no M32/M33/M35 process
holds the GPU.  Otherwise prints the reasons and exits 1."""
import os
import subprocess
import sys

O = os.path.expanduser('~/coding-rna/outputs')
FILES = [f'{O}/M35-PLANTCAD2-LORA/M35_DONE', f'{O}/M34-NTV2-SCALEUP/GO_STAGE1']     # GO_STAGE1 is created by the coordinator after the user's backbone decision; never by this code
PAT = ('M32-C1-PLANTCAD2', 'M33-NTV2-FOLLOWUP', 'M35-PLANTCAD2-LORA')
IGNORE_CMD = ('pgrep', 'probe_m32.sh', 'wait_m32.sh', 'm34_gate.py', 'm34_queue.sh', 'grep ', 'ps aux')


def mine():
    s, p = set(), os.getpid()
    while p > 1:
        s.add(p)
        try:
            p = int(open(f'/proc/{p}/stat').read().rsplit(')', 1)[1].split()[1])
        except Exception:
            break
    return s


def scan():
    bad, own = [], mine()
    for d in os.listdir('/proc'):
        if not d.isdigit() or int(d) in own:
            continue
        try:
            cmd = open(f'/proc/{d}/cmdline', 'rb').read().replace(b'\0', b' ').decode(errors='replace').strip()
            cwd = os.readlink(f'/proc/{d}/cwd')
        except Exception:
            continue
        if any(x in cmd for x in IGNORE_CMD) or cmd.startswith(('ssh ', 'sshd')):
            continue
        if any(x in cmd or x in cwd for x in PAT):
            bad.append((int(d), cmd[:160]))
    return bad


def gpu_pids():
    try:
        out = subprocess.run(['nvidia-smi', '--query-compute-apps=pid', '--format=csv,noheader'], capture_output=True, text=True, timeout=30).stdout
        return [int(x) for x in out.split() if x.strip().isdigit()]
    except Exception:
        return []


def main():
    reasons = [f'missing {f}' for f in FILES if not os.path.exists(f)]
    bad = scan()
    reasons += [f'process {p}: {c}' for p, c in bad]
    own_gpu = {p for p, _ in bad} & set(gpu_pids())
    reasons += [f'GPU process {p}' for p in own_gpu]
    if reasons:
        print('BLOCKED', ' | '.join(reasons))
        return 1
    print('OK')
    return 0


if __name__ == '__main__':
    sys.exit(main())
