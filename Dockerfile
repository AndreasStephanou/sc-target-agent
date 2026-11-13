# sc-target-agent: the frozen analysis environment plus the code, for reviewers and CI.
#
# The image is the environment. Run your own clone inside it so runs/ and data/ stay on your disk:
#   docker run --rm -v "$PWD":/work -w /work ghcr.io/andreasstephanou/sc-target-agent:latest python scripts/sctarget.py demo
#   apptainer exec --bind "$PWD":/work --pwd /work docker://ghcr.io/andreasstephanou/sc-target-agent:latest python scripts/sctarget.py demo
# Without a mount the copy baked into the image at /opt/sc-target-agent is used (tests, demo, check-env).
# The agent run needs the Claude Code CLI login or ANTHROPIC_API_KEY; neither is in the image.
#
# Package versions come from environment/*.lock (scripts/freeze_env.sh), not from environment.yml,
# so a build months from now installs the versions that produced the results tags.
FROM condaforge/miniforge3:latest

ENV PYTHONDONTWRITEBYTECODE=1 PYTHONUNBUFFERED=1 PIP_NO_CACHE_DIR=1 \
    MPLBACKEND=Agg NUMBA_CACHE_DIR=/tmp/numba_cache

# OS library the locked packages expect from the host: 8 of the 182 conda packages come from the
# defaults channel, _openmp_mutex among them, whose defaults build does not pull in a conda libgomp the
# way the conda-forge one does; without the system library `import sklearn` fails with
# "libgomp.so.1: cannot open shared object file".
RUN apt-get update && apt-get install -y --no-install-recommends libgomp1 && rm -rf /var/lib/apt/lists/*

# conda layer: explicit URLs + md5, installed without a solver
COPY environment/conda-linux-64.lock /tmp/conda-linux-64.lock
RUN conda create -y -p /opt/sctarget --file /tmp/conda-linux-64.lock && conda clean -afy

# pip layer: exact versions, no resolution
COPY environment/pip-linux-64.lock.txt /tmp/pip-linux-64.lock.txt
RUN /opt/sctarget/bin/pip install --no-deps -r /tmp/pip-linux-64.lock.txt

ENV PATH=/opt/sctarget/bin:$PATH

# the code
WORKDIR /opt/sc-target-agent
COPY . .
RUN pip install --no-deps -e . && mkdir -p runs cache

CMD ["python", "scripts/sctarget.py", "demo"]
