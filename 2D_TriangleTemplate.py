import turtle
import numpy as np
import pandas as pd
import random

nx = 10
ny = 10

u = np.ones((ny, nx))
placement = 8
elevation = 7
triangleStart = 0
base = 2
triangleHeight = elevation - base
u[base:triangleHeight, placement] = 0
u[base, triangleStart:placement] = 0
tilt = triangleHeight / (placement - triangleStart)
print(tilt)
for it in range(placement - triangleStart +1):
    u[base:int(base + it * tilt)+1, triangleStart + it] = 0

df = pd.DataFrame(u)
print(df)

turtlestart = -200
turtle.setup(1000, 550)
turtle.speed(0)  # Fastest speed

for m in range(ny):
    for b in range(nx):
        color = (0, 0,0)
        turtle.color(color)
        goto = turtlestart + b * 25
        turtle.penup()
        turtle.goto(goto, m * 25)
        turtle.pendown()



        # Draw the square
        for n in range(4):
            turtle.forward(25)
            turtle.right(90)

        turtle.penup()
        turtle.goto(turtlestart + b * 25 + 12.5, m * 25 - 20)
        if u[m, b] == 0:
            color = (1, 0, 0)
            turtle.color(color)
        turtle.write(str(f"{u[m, b]:.2f}"), align="center", font=("Arial", 8, "normal"))
    color = (0, 0, 1)
    turtle.color(color)
    turtle.penup()
    turtle.goto(turtlestart - 10, m * 25 - 20)
    turtle.write(str(m), align="center", font=("Arial", 8, "normal"))

for numberx in range(nx):
    gotox = turtlestart + numberx * 25
    turtle.penup()
    turtle.goto(gotox + 12.5, -45)
    turtle.write(str(numberx), align="center", font=("Arial", 8, "normal"))

# Hide the turtle
turtle.hideturtle()

# Keep the turtle window open until it is manually closed
turtle.done()
