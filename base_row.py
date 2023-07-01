import turtle
import numpy as np  # we're importing numpy


nx = 11
dx = 2 / (nx - 1)
nt = 50 # nt is the number of timesteps we want to calculate
dt = .0025  # dt is the amount of time each timestep covers (delta t)




# Set up the turtle screen
turtle.setup(1500, 200)
turtle.speed(0)  # Fastest speed


for i in range(nx):
    goto = -450 + i * 25
    turtle.penup()
    turtle.goto(goto, 0)
    turtle.pendown()

    # Draw the square
    for n in range(4):
        turtle.forward(25)
        turtle.right(90)

    turtle.penup()
    turtle.goto(goto+12.5, -45)
    turtle.write(str(i), align="center", font=("Arial", 8, "normal"))



# Hide the turtle
turtle.hideturtle()
# Keep the turtle window open until it is manually closed
turtle.done()
