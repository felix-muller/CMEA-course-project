from numpy import *
from pylab import *
import sys

basename = "square" #sys.argv[1]
errors = loadtxt("%s_errors.txt" % basename)
errorsH1 = loadtxt("%s_errorsH1.txt" % basename)
numberOfElements = loadtxt("%s_resolutions.txt" % basename)
print("L2:")
print(-diff(log(errors)) / diff (log(numberOfElements)))
print("EOC (L2): ", -polyfit(log(numberOfElements), log(errors), 1)[0])

print("H1:")
print(-diff(log(errorsH1)) / diff (log(numberOfElements)))
print("EOC (H1): ", -polyfit(log(numberOfElements), log(errorsH1), 1)[0])

plt.subplot(121)
loglog(numberOfElements, errors, "-o")
title("$L^2$ error")
xlabel("Number of free vertices")
ylabel("$||u_{FEM}-u||_{L^2}$")
grid("on")

plt.subplot(122)
title("$H^1$ error")
loglog(numberOfElements, errorsH1, "-o")

xlabel("Number of free vertices")
ylabel("$||u_{FEM}-u||_{H^1}$")
grid("on")

show()



# print(diff(log(errors)) / diff (log(numberOfElements)))
# print(polyfit(log(numberOfElements), log(errors), 1))
