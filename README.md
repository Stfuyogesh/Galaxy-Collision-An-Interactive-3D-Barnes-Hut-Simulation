# Galaxy Collision

### An Interactive 3D Galaxy Simulation

A browser-based simulation of two galaxies interacting through gravity. Built using HTML, CSS, and JavaScript, the project combines physics, algorithmic optimization, and real-time visualization.

## What is Barnes–Hut?

The Barnes–Hut algorithm is a method for calculating gravitational forces between many objects more efficiently.

Imagine simulating 4,000 stars. Calculating the gravitational influence between every pair requires millions of interactions in each force evaluation.

Barnes–Hut solves this problem by dividing 3D space into smaller regions using a structure called an **octree**. Instead of calculating the effect of every distant star individually, the algorithm can approximate a group of stars as one object located at their combined center of mass.

This reduces the computational workload while introducing some approximation error.

- **Brute force:** Calculates the gravitational influence of individual stars directly, requiring O(N²) work.
- **Barnes–Hut:** Groups distant stars to reduce unnecessary calculations, typically achieving approximately O(N log N) complexity for well-behaved distributions.

The simulation lets users adjust the accuracy and compare the computational workload of both approaches.

## How We Built It

**1. HTML: Structure**

We created the interface with canvas elements for the simulation and background, alongside controls for star count, settings, live statistics, and the introductory story.

**2. CSS: Visual Design**

We used CSS to create the dark space-themed interface, glowing visual elements, translucent panels, animations, and responsive layouts.

**3. JavaScript: Physics and Simulation**

JavaScript handles the entire simulation, including star positions, velocities, masses, gravitational forces, and camera movement.

We initialized two rotating stellar systems with different distributions and velocities, then calculated their gravitational interactions over successive simulation steps.

**4. Octree: Efficient Calculations**

We implemented a three-dimensional octree that recursively divides space into eight smaller regions. Each region stores information about the mass and center of mass of its contents.

The Barnes–Hut algorithm traverses this tree to determine when distant groups can be approximated and when individual regions need further inspection.

**5. Canvas: Rendering in 3D**

We projected the simulated 3D coordinates onto an HTML Canvas using a perspective camera. Stars are rendered with glowing sprites, while the galactic cores receive additional visual effects.

Users can rotate the camera, zoom in, follow individual stars, inspect performance metrics, and visualize the octree.

**6. Interactive Controls and Audio**

We added controls for star count, approximation accuracy, slow motion, and brute-force calculations. Ambient music is generated directly in the browser using the Web Audio API.

## Technology Stack

- HTML5
- CSS3
- Vanilla JavaScript
- HTML Canvas 2D
- Web Audio API
- Barnes–Hut algorithm
- 3D octree spatial partitioning

No frontend framework or build system is required.

## Limitations

This is an educational simulation, not a precision astrophysical model. It uses arbitrary simulation units, simplified initial conditions, gravitational softening, and approximate force calculations. Actual performance depends on the device and simulation settings.

---

**The goal:** demonstrate how an algorithmic optimization technique can make a computationally expensive gravitational simulation interactive in a web browser.