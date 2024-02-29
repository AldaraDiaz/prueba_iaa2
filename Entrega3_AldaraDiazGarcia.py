# RESOLUCIÓN DE 6.4.APLICACIÓN 2 --> ECUACIÓN DE ADVECCIÓN DE ONDA
# Aldara Díaz García

import math
from numpy import zeros
from scipy.linalg import lstsq
from matplotlib.pyplot import plot, show

import numpy as np 
import matplotlib.pyplot as plt

# Datos 
c=1 # Velocidad en m/s
L=3 
n=150 # Número de nodos
dx=L/n # Paso espacial
T=0.5 # Tiempo en segundos para la resolución

# MÉTODO EXPLÍCITO
maxk=1000 # Valor máximo de espacios en el tiempo para la resolución mediante el método explícito
dt=T/maxk # Paso temporal

# Defino las matrices
x=zeros((n+1,1)) # Vector que define la posición en el espacio
tiempo=zeros((maxk+1,1)) # Vector con todos los intantes de tiempo 
u=zeros((n+1,maxk+1)) # Matriz resultados

# Condiciones iniciales en t=0
for i in range(n+1):
    x[i]=i*dx
    u[i,0]=5+2*math.sin(math.pi*x[i,0])

# Condiciones de contorno para x=0
for k in range(maxk):
    u[0,k]=5
    tiempo[k]=k*dt

# Resolución empleando el método explícito
for k in range(maxk):
    for i in range(1,n):
        u[i,k+1]=u[i,k]-c*dt/(2*dx)*(u[i+1,k]-u[i-1,k])

print(c*dt/dx)

# Plotear solución empleando el método explícito
# Pasos de tiempo intermedios
t1=0.25*T
t2=0.5*T
t3=0.75*T

kt1=round(t1/dt)
kt2=round(t2/dt)
kt3=round(t3/dt)

plt.plot(x, u[:,0], "-o", x, u[:,kt1], "-s", x, u[:,kt2], "-+", x, u[:,kt3], "-*",x, u[:,maxk-1], "-h")
plt.xlabel('x')
plt.ylabel('u')
plt.title('Método explícito')
plt.show()

# MÉTODO DE LAX
maxk_lax=1000 # Valor máximo de espacios en el tiempo para la resolución mediante el esquema de Lax
dt=T/maxk_lax # Paso temporal

# Defino las matrices
x=zeros((n+1,1)) # Vector que define la posición en el espacio
tiempo=zeros((maxk_lax+1,1)) # Vector con todos los intantes de tiempo 
u=zeros((n+1,maxk_lax+1)) # Matriz resultados

# Condiciones iniciales en t=0
for i in range(n+1):
    x[i]=i*dx
    u[i,0]=5+2*math.sin(math.pi*x[i,0])

# Condiciones de contorno para x=0
for k in range(maxk_lax):
    u[0,k]=5
    tiempo[k]=k*dt

# Resolución empleando el esquema de Lax
for k in range(maxk_lax):
    for i in range(1,n):
        u[i,k+1]=1/2*(u[i+1,k]+u[i-1,k])-c*dt/(2*dx)*(u[i+1,k]-u[i-1,k])

# Plotear solución empleando el esquema de Lax
# Pasos de tiempo intermedios
t1=0.25*T
t2=0.5*T
t3=0.75*T

kt1=round(t1/dt)
kt2=round(t2/dt)
kt3=round(t3/dt)

plt.plot(x, u[:,0], "-o", x, u[:,kt1], "-s", x, u[:,kt2], "-+", x, u[:,kt3], "-*", x, u[:,maxk_lax-1], "-h")
plt.xlabel('x')
plt.ylabel('u')
plt.title('Método de Lax')
plt.show()