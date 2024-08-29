# Finite-volume-method-for-solving-Laplace-equation-on-a-distorted-annular-domain-with-Cartesian-grid
Projects completed as part of the Advanced Computational Fluid Dynamics (ME634) course at IIT Kanpur. The project encompass various simulations and analyses, demonstrating proficiency in numerical methods, turbulence modeling, and computational techniques in fluid dynamics.


This project demonstrates the application of the Finite Volume Method (FVM) to solve Laplace’s equation
 ∂2ϕ∂x2+∂2ϕ∂y2=0\frac{\partial^2 \phi}{\partial x^2} + \frac{\partial^2 \phi}{\partial y^2} = 0∂x2∂2ϕ+∂y2∂2ϕ=0 
on an annular domain between two concentric circles. The problem is solved using Cartesian coordinates, showcasing the FVM’s ability to handle distorted domains while benefiting from a known exact solution.
The provided matlab code is used to perform the numerical simulations.
Computational Domain
The computational domain is an annular region defined by the following parameters:
•	Inner Radius rW=0.1r_W = 0.1rW=0.1
•	Outer Radius rZ=1.0r_Z = 1.0rZ=1.0
•	Angular Start θX=0∘\theta_X = 0^\circθX=0∘
•	Angular End θY=90∘\theta_Y = 90^\circθY=90∘
•	Over-relaxation Parameter λ=1.5\lambda = 1.5λ=1.5
Boundary conditions for Laplace’s equation are:
•	ϕ=0\phi = 0ϕ=0 on the bottom horizontal boundary WX
•	ϕ=rsin⁡(θ)\phi = r \sin(\theta)ϕ=rsin(θ) on the left vertical boundary XY
•	ϕ=1\phi = 1ϕ=1 on the top horizontal boundary YZ
•	ϕ=rsin⁡(θ)\phi = r \sin(\theta)ϕ=rsin(θ) on the right vertical boundary WZ
The exact solution for validation is ϕ=rsin⁡(θ)\phi = r \sin(\theta)ϕ=rsin(θ).
