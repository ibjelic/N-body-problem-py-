# N-body problem python
2D N body problem in python with electric and gravitational force

Graphical preview is not the best but it works. 
Class body contains all info about body and interactions with other bodies, and as arguments most of functions take body[] (array).

System for preventing detection of collision works in a way that calculates speeds after it and moves bodies with that speed untill there is no collision, so second body wont detect collision and change parameters again.

Simulation function gives back a matrix of data (info * bodies * step, and info is [x,y,Vx,Vy]), after simulation matplotlib animates data and shows animation.

Screenshot of live graph:
![Screenshot of live graph](https://imgur.com/i7EqZIU.png)

Graphs are animations of positions and speeds over time (time of simulation)


