---
name: graphify
description: Use Graphify to analyze code dependencies, call chains, file references, and project structure before modifying MATLAB bearing simulation code. Use this skill when the user asks to understand, trace, inspect, or modify the repository safely.
---

# Graphify skill for MATLAB bearing simulation

Use this skill to analyze the repository structure before making code changes.

## Purpose

Graphify is used only for:
- dependency analysis
- call-chain analysis
- file-reference tracing
- code-structure understanding
- identifying risky modification points

It must not be used to validate physical formulas or replace MATLAB smoke tests.

## Hard rules

1. Do not modify physical formulas unless explicitly requested.
2. Do not change contact mechanics, oil-film, stiffness, dynamics, load, impurity, or bearing geometry models unless explicitly requested.
3. Do not scan temp_workspace, outputs, *.mat, *.log, *.fig, *.png, *.sqlite, or other generated files.
4. Do not generate or open graph.html unless explicitly requested.
5. Prefer text-only dependency summaries.
6. Do not run MATLAB during Graphify-only analysis.
7. Do not use historical outputs as new results.

## Required analysis targets

When analyzing this repository, focus on:

1. run_all_cases_with_reports call chain.
2. temp_workspace and worker01 references.
3. result111.mat, result222.mat, result333.mat, result444.mat read/write chain.
4. zz1.mat, zz2.mat, Jzz.mat, Q1nonload.mat, Ph2ii.mat read/write chain.
5. completed=true assignment locations.
6. failure_report generation locations.
7. www1, loadi, loadj data flow in the roller bearing workflow.
8. ball bearing and roller bearing branch points.

## Output format

For dependency analysis, output only:

1. Main call chain.
2. MAT-file communication chain.
3. Concurrency risk points.
4. Roller state synchronization risk points.
5. completed=true misclassification risk points.
6. Recommended modification files.

Do not output images, graph visualizations, or long logs.
