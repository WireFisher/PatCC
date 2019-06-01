#!/bin/bash

kill -9 `ps -ef | grep run_test | awk '{print $2}'`
kill -9 `ps -ef | grep run_all_test | awk '{print $2}'`
