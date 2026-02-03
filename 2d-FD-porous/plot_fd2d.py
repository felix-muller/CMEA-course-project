import mpl_toolkits.mplot3d
from numpy import *
from pylab import *

u = loadtxt("u_fd.txt")
N = int(sqrt(u.shape[0]))
print("N=%d" % N)
u.shape = (N,N)

fig = figure()
ax = fig.add_subplot(111, projection='3d')
x,y = mgrid[0:1:N*1j, 0:1:N*1j]
ax.plot_surface(x,y,u)
xlabel('x')
ylabel('y')
#zlabel('u')
title('Approximation to the solution')
show()

contourf(x,y,abs(u - sin(2*pi*x)*sin(2*pi*y)))
colorbar()
xlabel('x')
ylabel('y')
title('Approximation error')
show()

plot(x[:, 50], u[:, 50], '-*', label='$U_{i,50}$')
plot(x[:,50], sin(2*pi*y[:,50])*sin(2*pi*x[:,50]), '--',label='exact')
legend()
show()
