#!/usr/bin/env bash
cd $(dirname "$0")
export PYTHONPATH=src
python tests/test_time.py
