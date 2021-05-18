FROM python:3.8.3

WORKDIR /build/gurobi_installer

RUN wget -q https://packages.gurobi.com/9.1/gurobi9.1.0_linux64.tar.gz && \
    tar -xf gurobi9.1.0_linux64.tar.gz -C /usr/share

WORKDIR /build/switch_installer

RUN git clone https://github.com/switch-model/switch.git
WORKDIR /build/switch_installer/switch
RUN    pip install --upgrade --editable .

COPY examples/wecc/ ./wecc/

ENV PATH="$PATH:/usr/share/gurobi910/linux64/bin" \
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:/usr/share/gurobi910/linux64/lib \
    GUROBI_HOME='/usr/share/gurobi910/linux64' \
    GRB_LICENSE_FILE='/usr/share/gurobi_license/gurobi.lic' \
