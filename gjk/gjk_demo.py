import pygame
import numpy as np
from gjk_algo import GJK, Circle


def Cube():
    glBegin(GL_LINES)
    for edge in edges:
        for vertex in edge:
            glVertex3fv(vertices[vertex])
    glEnd()

pygame.init()

screen_width = 640
screen_height = 480
grid_size = 10.0

size = [screen_width, screen_height]
screen = pygame.display.set_mode(size)
pygame.display.set_caption('GJK Algorithm in 2D')

done = False
clock = pygame.time.Clock()

WHITE = (255, 255, 255)
BLACK = (0, 0, 0)
RED = (255, 0, 0)
ORANGE = (200, 128, 0)
GRAY = (64, 64, 0)
PURPLE = (128, 128, 0)
GREEN = (0, 255, 0)
BLUE = (0, 0, 255)


def s2c(x, y=None):
    if y is None: x, y = x
    origin_x = screen_width * 0.5
    origin_y = screen_height * 0.5
    cx = (x - origin_x) / grid_size
    cy = (origin_x - y) / grid_size
    return (cx, cy)


def c2s(x, y=None):
    if y is None: x, y = x
    origin_x = screen_width * 0.5
    origin_y = screen_height * 0.5
    sx = origin_x + x * grid_size
    sy = origin_y - y * grid_size
    return (int(sx), int(sy))


polygon1 = Circle((12, 12), 2) #[(1,4), (5,3), (10, 5), (11,9), (5,9), (2,8)]
polygon2 = [(10, 10), (12, 6), (16, 8), (15, 11)]

## don't calculate minkowsky distance of a circle
# minkowsky_diff = []
# for p1 in polygon1:
#     for p2 in polygon2:
#         minkowsky_diff.append((p1[0]-p2[0], p1[1]-p2[1]))

algo = GJK(polygon1, polygon2, init_d=(-1, 0))
status = algo.step()
print "init step:", status

def drawLines(screen, color, close, c_shape):
    pygame.draw.lines(screen, color, True, [c2s(p) for p in c_shape])

def drawShape(screen, color, shape):
    if isinstance(shape, Circle):
        pygame.draw.circle(screen, color, c2s(shape.center), int(shape.radius * grid_size), 1)
    else:
        drawLines(screen, color, True, shape)

while not done:
    clock.tick(24)

    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            done = True
        if event.type == pygame.KEYDOWN:
            if event.key == pygame.K_SPACE:
                if status == "In progress":
                    status = algo.step()
                    print status
                else:
                    done = True

    screen.fill(WHITE)

    origin_rect = [screen_width/2-2, screen_height/2-2, 4, 4]
    pygame.draw.rect(screen, BLACK, origin_rect)
    drawShape(screen, PURPLE, polygon1)
    drawShape(screen, BLUE, polygon2)
    # for p in minkowsky_diff:
    #     pygame.draw.circle(screen, GREEN, c2s(p), 2)


    # draw old simplex
    if len(algo.simplex) > 3:
        old_simplex = algo.simplex[:-1]
        drawLines(screen, ORANGE, True, old_simplex)
    
        new_simplex = algo.simplex[-3:]
        drawLines(screen, RED, True, new_simplex)
    else:
        drawLines(screen, RED, False, algo.simplex)


    nextNormal = algo.d
    pygame.draw.line(screen, GRAY, c2s(0, 0), c2s(nextNormal))
    
    pygame.display.flip()

pygame.quit()
