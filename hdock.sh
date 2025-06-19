#!/bin/bash

# 开始执行脚本
echo "Running PPbench affinity evaluation..."
nohup python /root/autodl-tmp/PP_generate_v1/utils/evaluate/hdock_PPbench.py > /root/autodl-tmp/PP_generate_v1/data/output/PPbench_vina.log 2>&1 &
echo "PPbench evaluation done. Waiting for one minute..."
sleep 60

echo "Running RFdiffusion affinity evaluation..."
nohup python /root/autodl-tmp/PP_generate_v1/utils/evaluate/hdock_RFdiffusion.py > /root/autodl-tmp/PP_generate_v1/data/output/RFdiffusion_vina.log 2>&1
echo "RFdiffusion evaluation done. Waiting for one minute..."
sleep 60

echo "Running ppflow affinity evaluation..."
nohup python /root/autodl-tmp/PP_generate_v1/utils/evaluate/hdock_PPFlow.py > /root/autodl-tmp/PP_generate_v1/data/output/ppflow_vina.log 2>&1
echo "ppflow evaluation done. Waiting for one minute..."
sleep 60

echo "Running peptide generation affinity evaluation..."
nohup python /root/autodl-tmp/PP_generate_v1/utils/evaluate/hdock.py > /root/autodl-tmp/PP_generate_v1/data/output/generate_vina.log 2>&1
echo "Peptide generation evaluation done."

echo "All evaluations completed successfully."
