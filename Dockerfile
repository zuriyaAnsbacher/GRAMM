FROM python:3.8

ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

ARG APP_DIR="/app"

WORKDIR $APP_DIR

RUN python -m pip install --upgrade pip 

COPY app/requirements.txt .

RUN pip install -r requirements.txt

COPY app .

# RUN conda install -c anaconda mkl-service
# RUN conda install -c conda-forge mkl_fft
# RUN conda install -c anaconda mkl_random

# RUN pip install example --use-feature=2020-resolver

CMD python server.py
