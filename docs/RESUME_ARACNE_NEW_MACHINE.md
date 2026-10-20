# Resume the v2 ARACNe build on a new machine (MAIL-15)

Written 2026-10-20. The v2 network build (`scripts/43_run_aracne.sh`) was paused on MAIL-10 (128.2.112.20) on
2026-10-19 and moves to MAIL-15 (128.2.112.22, Lambda #6). See `docs/PROGRESS.md` (2026-10-19) for the pause.

**State when paused:** COADREAD finished (323,571 edges); BRCA 68 / 100 and BRCA_LUM 84 / 100 bootstraps; the other
26 networks and the TF-only set (`REGSET=tf`) not started. The script resumes from the last finished bootstrap and
skips finished networks. Output depends on the seed, not the thread count, so THREADS / PARALLEL change speed only.

## 1. SSH config (on your laptop)

Add to `~/.ssh/config`:

```
# MAIL-15
# Lambda #6 (3 x Ada A6000)
Host 15
    HostName 128.2.112.22
    ForwardAgent yes
    User kirti
```

Then:

```bash
ssh 15            # first time: accept the host key, log in with your password
ssh-copy-id 15    # optional: install your laptop key so no password is needed
```

Off campus you may need the CMU VPN first.

## 2. Check the machine (on MAIL-15)

The GPUs do not matter: ARACNe-AP is CPU-only. Cores and RAM decide THREADS and PARALLEL (each network in parallel
needs up to ~60 GB heap, `-Xmx60G` in the script).

```bash
nproc; lscpu | grep -E "Model name|Core|Socket"; free -g; df -h ~
```

## 3. Get the code (on MAIL-15)

`ForwardAgent yes` passes your laptop's GitHub key through:

```bash
ssh -T git@github.com    # should greet you by username
mkdir -p ~/projects && cd ~/projects
git clone git@github.com:kirtihargunani-cmu/repurposing.git
```

If that fails: `gh auth login`, then clone over HTTPS.

## 4. Conda, if missing (on MAIL-15)

```bash
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p ~/miniforge3
```

## 5. Java + ARACNe-AP, same versions as before (on MAIL-15)

OpenJDK 17.0.18, ARACNe-AP commit 30967ea:

```bash
cd ~/projects/repurposing && bash scripts/07_setup_aracne.sh
```

## 6. Copy the inputs and finished bootstraps (~2.6 GB, gitignored)

On MAIL-15:

```bash
mkdir -p ~/projects/repurposing/data/processed/step5 ~/projects/repurposing/logs
```

On MAIL-10. MAIL-10 can reach MAIL-15 but has no key there: either log in to MAIL-10 with `ssh -A` (your laptop key,
installed on MAIL-15 in step 1, is forwarded), or add MAIL-10's `~/.ssh/id_*.pub` to `~/.ssh/authorized_keys` on
MAIL-15.

```bash
cd ~/projects/repurposing
rsync -avP data/processed/step5/aracne_inputs_v2 data/processed/step5/aracne_v2 \
  kirti@128.2.112.22:~/projects/repurposing/data/processed/step5/
```

Check on MAIL-15: `ls data/processed/step5/aracne_v2/BRCA | grep -c '^bootstrapNetwork'` should print 68
(BRCA_LUM 84; COADREAD has `network.txt`).

## 7. Resume the build (on MAIL-15)

Set THREADS (per network) and PARALLEL (networks at once) from step 2, e.g. 16 and 1 for 16 threads; keep
PARALLEL x 60 GB under the RAM. Full set first, then the TF-only set:

```bash
cd ~/projects/repurposing
nohup setsid bash -c 'THREADS=16 PARALLEL=1 VERSION=v2 scripts/43_run_aracne.sh >> logs/aracne_v2.log 2>&1; THREADS=16 PARALLEL=1 VERSION=v2 REGSET=tf scripts/43_run_aracne.sh > logs/aracne_v2_tf.log 2>&1' > /dev/null 2>&1 < /dev/null &
```

Progress: `tail logs/aracne_v2.log`; bootstraps per network:
`for d in data/processed/step5/aracne_v2/*/; do echo "$(basename $d) $(ls $d | grep -c '^bootstrapNetwork')"; done`.

To pause safely, kill the whole process group (not only java, or the loop moves on to the next network):
`kill -TERM -- -<pgid>` (find it with `ps -eo pid,pgid,cmd | grep 43_run_aracne`).

## After the build

Copy `data/processed/step5/aracne_v2*` back (or continue there), then: v2 regulons with cross-mappable edges removed,
VIPER v2, checks C2, C5, C7-C10, C12, and the BRCA_LUM vs BRCA comparison (D68; fix the "material change" cutoff
first). See `docs/SPEC_rebuild_v2.md`.
