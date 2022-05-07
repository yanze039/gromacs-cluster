FROM dp-harbor-registry.cn-zhangjiakou.cr.aliyuncs.com/dp/hermite-fold-cpu:1.0.base
SHELL ["/bin/bash", "-c"]
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys A4B469963BF863CC
RUN apt-get update
RUN echo -e "6\n1\n" | apt-get install -y gcc g++ cmake
ENV CC=/usr/bin/gcc
ENV CXX=/usr/bin/g++

ADD cmake-3.16.0-Linux-x86_64.tar.gz /root/
RUN cp -R /root/cmake-3.16.0-Linux-x86_64/bin /usr/ \
    && cp -R /root/cmake-3.16.0-Linux-x86_64/share /usr/

# RUN wget -O root/gromacs-2021.4.tar.gz https://github.com/gromacs/gromacs/archive/refs/tags/v2021.4.tar.gz --no-check-certificate
ENV GMX_PREFIX "/opt/gromacs-2021.4"
ADD gromacs-2021.4.tar.gz /root/
RUN cd /root/gromacs-2021.4 \
    && mkdir build \
    && cd build \
    && cmake .. -DGMX_BUILD_OWN_FFTW=ON \
                -DCMAKE_INSTALL_PREFIX=$GMX_PREFIX \
                -DGMX_GPU=OFF  \
                -DGMX_SIMD=avx_512 \
    && make -j 8 \
    && make install
RUN rm -rf /root/gromacs-2021.4 /root/gromacs-2021.4.tar.gz

# environment variables set manually since github actions will overwrite entrypoints
ENV GROMACS_DIR=$GMX_PREFIX
ENV GMXBIN=$GROMACS_DIR/bin
ENV GMXDATA=$GROMACS_DIR/share/gromacs
ENV GMXLDLIB=$GROMACS_DIR/lib
ENV GMXMAN=$GROMACS_DIR/share/man
ENV LD_LIBRARY_PATH=$GMXLDLIB:$LD_LIBRARY_PATH
RUN echo "source $GMX_PREFIX/bin/GMXRC" >> /root/.bashrc
SHELL ["/bin/bash", "-c"]
# ADD test.pdb /root/
# RUN source $GMX_PREFIX/bin/GMXRC \
#     && which gmx \
#     && cd /root/ \
#     && echo "3\n" | gmx cluster -f test.pdb -s test.pdb -method gromos -cutoff 0.15

