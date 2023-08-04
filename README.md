# Title of your Project

This repository contains code and related files for [your project name]. The project involves utilizing TLE (Two-Line Element Set) data for various applications, including orbit prediction, visibility check, and conflict resolution.

## Table of Contents
1. [Introduction to TLE](#introduction-to-tle)
2. [Problem Setup](#problem-setup)
3. [Rough Propagation](#rough-propagation)
4. [Visibility Check](#visibility-check)
5. [Conflict Resolution](#conflict-resolution)
6. [Precise Propagation](#precise-propagation)
7. [Conclusion](#conclusion)

## Introduction to TLE

TLE, or Two-Line Element set, is a data format used by NORAD and NASA to represent an Earth Satellite's orbit. It contains information such as the inclination, right ascension, and mean anomaly which are crucial for tracking a satellite's movement and position in space. Our project leverages this data for various operations related to satellite management.

## Problem Setup

The primary challenge our project addresses is the management of satellite operations. This involves predicting the satellite's orbit, ensuring its visibility from a specific ground station at a given time, and resolving any potential conflicts with other space bodies. All of these tasks heavily depend on the TLE data of the satellite.

## Rough Propagation

In our project, we first perform a rough or coarse orbit prediction. This is done using a simplified model to estimate the satellite's position and velocity, which, although less accurate, is less computationally intensive. This initial prediction serves as a stepping stone for the upcoming operations.

## Visibility Check

Next, we use the TLE data to perform a visibility check. This process verifies when our satellite is visible or 'in view' from a specified ground station. It is crucial for planning communication or observational operations with the satellite.

## Conflict Resolution

One of the primary concerns in satellite management is conflict resolution. As space becomes more crowded, it is crucial to predict potential conjunctions or close approaches with other objects. By leveraging TLE data, we can plan manoeuvres to avoid possible collisions, thus maintaining the satellite's safety.

## Precise Propagation

When a higher degree of precision is required, we perform precise propagation. Unlike rough propagation, this involves the use of more accurate physics models and factors in elements overlooked in the initial prediction. It is computationally more intensive but offers a more accurate prediction, which is essential for mission-critical tasks.

## Conclusion

Our project utilizes TLE data to navigate the challenges of satellite operations effectively. From predicting orbits and checking visibility to resolving conflicts, we have developed a robust and efficient system that ensures smooth and safe satellite operations.
