# GNC Simulation Environment

This is a GNC simulation environment that lets you simulate flight software running on a satellite, complete with attitude dynamics, orbit dynamics, hardware models, etc.

### Table of Contents  
 &nbsp;&nbsp;&nbsp;&nbsp; [Features](#Features)  

 &nbsp;&nbsp;&nbsp;&nbsp; [How to Use](#How-to-Use)  
 &nbsp;&nbsp;&nbsp;&nbsp; [Setting up Physics Engine](#Setting-up-Physics-Engine)  
 &nbsp;&nbsp;&nbsp;&nbsp; [Setting up Flight Software](#Setting-up-Flight-Software)  

 &nbsp;&nbsp;&nbsp;&nbsp; [Running the Simulation](#Running-the-Simulation)  
 &nbsp;&nbsp;&nbsp;&nbsp; [Running Example Simulation - Attitude Dynamics](#Running-the-Simulation-Attitude-Example)  
 &nbsp;&nbsp;&nbsp;&nbsp; [Running Example Simulation - Orbit Dynamics](#Running-the-Simulation-Orbit-Example)  


## Features

### Rigid Body Simulation 

  Attitude is determined using quaternions, propagated according to Euler's equation of rigid bodies. Working on a vibe model but isn't a priority at the moment.

### Orbit Propagator

  Orbit propagator handles classic orbit elements, r and v vectors, and accepts control and disturbance forces. Has a built-in drag model but you can use whatever model you want, stick it in the propagator as a disturbance force.

### GN&C Functions

  Pre-built GN&C functions such as Kalman filter, sliding mode controller, feedforward controller, pid controller, and more. Also implemented the Q-Law guidance equations for orbit transfers. You can make your own and simulate how they operate in the physics sim.

### Hardware Models

  

### Visualizer

  

## How to Use

## Setting up Physics Engine

## Setting up Flight Software

## Running the Simulation

## Running Example Simulation - Attitude Dynamics

## Running Example Simulation - Orbit Dynamics
