import pygame
from pygame.locals import *

import OpenGL
from OpenGL.GL import *
from OpenGL.GLU import *

import numpy as np
from gjk_algo import GJK3, getMinkowskyDiff

vertices= (
    (1, -1, -1),
    (1, 1, -1),
    (-1, 1, -1),
    (-1, -1, -1),
    (1, -1, 1),
    (1, 1, 1),
    (-1, -1, 1),
    (-1, 1, 1)
)

edges = (
    (0,1),
    (0,3),
    (0,4),
    (2,1),
    (2,3),
    (2,7),
    (6,3),
    (6,4),
    (6,7),
    (5,1),
    (5,4),
    (5,7)
)


def Cube(vertices, color=(0, 0, 0)):
    glColor3f(color[0] * 1.0 / 255, color[1] * 1.0 / 255, color[2] * 1.0 / 255)
    glBegin(GL_LINES)
    for edge in edges:
        for vertex in edge:
            glVertex3fv(vertices[vertex])
    glEnd()

def drawVertexCloud(verts):
    glPointSize(2.0)
    glColor3f(0, 1, 0)
    glBegin(GL_POINTS)
    for v in verts:
        glVertex3f(v[0], v[1], v[2])
    glEnd()

def drawTetrahedron(verts):
    glColor3f(1, 0, 0)
    glBegin(GL_LINES)
    for i in range(len(verts)):
        for j in range(i+1, len(verts)):
            glVertex3fv(verts[i])
            glVertex3fv(verts[j])
    glEnd()

def runProgram(func_list):
    for f in func_list:
        f()

def v4(v3, i):
    return np.array([v3[0], v3[1], v3[2], i])


done = False

WHITE = (255, 255, 255)
BLACK = (0, 0, 0)
RED = (255, 0, 0)
ORANGE = (200, 128, 0)
GRAY = (64, 64, 0)
PURPLE = (128, 0, 128)
GREEN = (0, 255, 0)
BLUE = (0, 0, 255)

KEY_PROCEED = pygame.K_SPACE
KEY_EDIT_MODE = pygame.K_1
KEY_VIEW_MODE = pygame.K_2
KEY_LEFT = pygame.K_a
KEY_RIGHT = pygame.K_d

EDIT_MODE = 1
VIEW_MODE = 2

x_rot = 0.0
y_rot = 0.0

translate_1 = np.random.rand(3) * 6.0 - 2.5
translate_2 = np.random.rand(3) * 6.0 - 3.5


# pygame settings
screen_width = 640
screen_height = 480
grid_size = 10.0

size = [screen_width, screen_height]
screen = pygame.display.set_mode(size, DOUBLEBUF|OPENGL)
pygame.display.set_caption("GJK Algorithm in 3D")

pygame.init()
clock = pygame.time.Clock()

# get some data cheap
glMatrixMode(GL_MODELVIEW)
glLoadIdentity()
glPushMatrix()
glRotatef(12, 1, 2, 3)
t_ = np.diagflat([1]*4)
t_[:3, 3] = translate_1.T
mvp1 = np.array(glGetFloatv(GL_MODELVIEW_MATRIX)).astype(np.float)
mvp1 = np.matmul(t_, mvp1)
glPopMatrix()

glPushMatrix()
glRotatef(5, -1, 0, 3)
t_ = np.diagflat([1]*4)
t_[:3, 3] = translate_2.T
mvp2 = np.array(glGetFloatv(GL_MODELVIEW_MATRIX)).astype(np.float)
mvp2 = np.matmul(t_, mvp2)
glPopMatrix()

# static gl settings
glMatrixMode(GL_PROJECTION)
glLoadIdentity()
gluPerspective(60, (1.0*size[0]/size[1]), 0.1, 50.0)
glTranslatef(0, 0, -10)
glDisable(GL_CULL_FACE)

glMatrixMode(GL_MODELVIEW)

cube1_t = [np.matmul(mvp1, v4(v, 1)) for v in vertices]
cube2_t = [np.matmul(mvp2, v4(v, 1)) for v in vertices]

algo = GJK3([v[:3] for v in cube1_t], [v[:3] for v in cube2_t])
status = algo.step()
print "Initial step:", status

# game loop
while not done:
    clock.tick(24)

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            done = True
        if event.type == pygame.KEYDOWN:
            if event.key == KEY_PROCEED:
                if status == "In progress":
                    status = algo.step()
                    print status
                else:
                    done = True
            if event.key == pygame.K_a:
                y_rot -= 15
            if event.key == pygame.K_d:
                y_rot += 15
            if event.key == pygame.K_w:
                x_rot -= 15
            if event.key == pygame.K_s:
                x_rot += 15

    glClearColor(1, 1, 1, 1)
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

    glPushMatrix()
    glRotatef(x_rot, 1, 0, 0)
    glRotatef(y_rot, 0, 1, 0)

    Cube([v[:3] for v in cube1_t], ORANGE)
    Cube([v[:3] for v in cube2_t], PURPLE)
    # drawVertexCloud(cube1_t)
    # drawVertexCloud(cube2_t)
    drawVertexCloud(getMinkowskyDiff(cube1_t, cube2_t))
    drawTetrahedron(algo.simplex)
    glPopMatrix()

    pygame.display.flip()


pygame.quit()
