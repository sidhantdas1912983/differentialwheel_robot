The desired trajectory was found out ussing RRT* but the file was not included in this as its already included in the 4Dof manipulator repro.
Once the trajectory was found out using navigation algorithm it was manually fed into the desired trajectory of the controller
one can change the desired trajectory as per the trajectory they desire and can visualise the movement using the animation algorithm present towards the end of the code.
Howevwer with changed trajectories the mean error of the controller might increase for which the Kp,KI and Kd parameter may be required to be tuned again.
