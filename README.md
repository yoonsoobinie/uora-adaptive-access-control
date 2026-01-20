# uora-adaptive-access-control

MATLAB simulator for adaptive access control in IEEE 802.11ax UORA using per-station transmission history.

## Overview
MATLAB simulation code for **adaptive access control** in IEEE 802.11ax UORA,
based on **per-station transmission history** (local observations only).

**Paper title:** Adaptive Access Control for IEEE 802.11ax UORA Based on Per-Station Transmission History  
**Authors:** Soobin Yoon, Eun-Chan Park (Dongguk University)

## What this repo contains
- A UORA sweep simulator that measures:
  - Throughput (Mbps)
  - Jain’s fairness index (success-count based)
  - Attempt rate
  - Average idle RA-RUs per TF
- The simulator follows standard UORA access procedure (OBO decrement per TF, uniform RU selection on attempt),
  while adaptively updating:
  - access threshold parameter (alpha, α)
  - OFDMA contention window (OCW)

## Key idea (high level)
Standard BEB-based OCW update can be insensitive to fast congestion changes, leading to either:
- more collisions, or
- more idle RA-RUs (under-utilization)

This approach estimates (per station) collision probability and waiting/idle probability from a sliding window,
then applies **sigmoid-normalized sensitivity** to adaptively update α and OCW without using global statistics.

## Reported results (from the paper)
- Average throughput improvement vs standard UORA: **~15.1%**
- Maximum throughput improvement:
  - (OCWmin, OCWmax) = (31, 511): **~50.0%**
  - (OCWmin, OCWmax) = (63, 1023): **~56.8%**
- Average idle RA-RUs decreased by about **~20%** (both settings in the reported experiments)

(Exact numbers depend on the simulation configuration.)

## Requirements
- MATLAB **R2022b or later** (uses local functions at the end of a script)

## How to run
1) Open MATLAB in the repository folder  
2) Run: `run_uora_adaptive_access.m`
