FROM ubuntu:24.04

# 이미지 제작 시 사용할 설정
ARG BUILD_JOBS=2
ARG SEAL_COMMIT=a4a428c0d630a31a2e4179fd964019f2e3b3264d

# 운영체제 패키지 설치
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        cmake \
        git \
        ca-certificates \
        python3 \
        python3-dev \
        python3-venv \
    && rm -rf /var/lib/apt/lists/*

# 컨테이너 전용 Python 가상환경
RUN python3 -m venv /opt/venv

ENV PATH="/opt/venv/bin:${PATH}" \
    CMAKE_PREFIX_PATH="/opt/seal" \
    PYTHONPATH="/app/build:/app/python" \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/matplotlib \
    OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1

# 검증한 SEAL 커밋을 내려받아 빌드·설치
RUN git clone https://github.com/CPET-lab/CPET_SEAL_Light.git \
        /opt/CPET_SEAL_Light \
    && git -C /opt/CPET_SEAL_Light checkout --detach "${SEAL_COMMIT}" \
    && cmake \
        -S /opt/CPET_SEAL_Light \
        -B /opt/seal-build \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX=/opt/seal \
        -DSEAL_USE_GAUSSIAN_NOISE=OFF \
        -DSEAL_USE_INTEL_HEXL=OFF \
        -DSEAL_BUILD_EXAMPLES=OFF \
        -DSEAL_BUILD_TESTS=OFF \
        -DSEAL_BUILD_BENCH=OFF \
    && cmake --build /opt/seal-build --parallel "${BUILD_JOBS}" \
    && cmake --install /opt/seal-build

# 이후 명령의 기준 디렉터리
WORKDIR /app

# Python 의존성 설치
COPY requirements.txt ./
RUN python -m pip install --no-cache-dir -r requirements.txt

# 배포 코드와 입력 캐시 복사
COPY CMakeLists.txt ./
COPY cpp/ ./cpp/
COPY python/ ./python/
COPY test/ ./test/
COPY cache/ ./cache/

# 배포 코드 빌드
RUN cmake -S . -B build \
        -DCMAKE_BUILD_TYPE=Release \
        -DPython_EXECUTABLE=/opt/venv/bin/python \
        -Dpybind11_DIR="$(python -m pybind11 --cmakedir)" \
    && cmake --build build \
        --target EREExperiment EBModule \
        --parallel "${BUILD_JOBS}"

# 프로젝트와 포함된 외부 코드의 라이선스 고지
COPY LICENSE THIRD_PARTY_NOTICES ./

# 명령을 지정하지 않았을 때의 기본 동작
CMD ["bash"]
