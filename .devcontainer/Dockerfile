FROM python:3.11-slim

ARG USERNAME=bafqmc
ARG USER_UID=1000
ARG USER_GID=1000

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bash \
        build-essential \
        ca-certificates \
        git \
        gfortran \
        openmpi-bin \
        libopenmpi-dev \
        libblas-dev \
        liblapack-dev \
    && rm -rf /var/lib/apt/lists/* \
    && groupadd --gid "${USER_GID}" "${USERNAME}" \
    && useradd --uid "${USER_UID}" --gid "${USER_GID}" \
        --create-home --shell /bin/bash "${USERNAME}" \
    && mkdir -p /workspaces/BAFQMC \
    && chown "${USERNAME}:${USERNAME}" /workspaces/BAFQMC

COPY benchmarks/paper/requirements.txt benchmarks/paper/requirements-ed.txt /tmp/bafqmc-requirements/

USER ${USERNAME}
ENV VIRTUAL_ENV=/home/${USERNAME}/.venv \
    PATH=/home/${USERNAME}/.venv/bin:$PATH \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PIP_NO_CACHE_DIR=1 \
    MPLBACKEND=Agg \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    NUMBA_NUM_THREADS=1

RUN python -m venv "${VIRTUAL_ENV}" \
    && python -m pip install \
        -r /tmp/bafqmc-requirements/requirements.txt \
        -r /tmp/bafqmc-requirements/requirements-ed.txt \
    && python -m pip check

WORKDIR /workspaces/BAFQMC
CMD ["bash"]
