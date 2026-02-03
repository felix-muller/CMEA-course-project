basename = 'square';

errors = load([basename, '_errors.txt']);
errorsH1 = load([basename, '_errorsH1.txt']);
numberOfElements = load([basename, '_resolutions.txt']);
disp("L2")
diffs = diff(-log(errors)) ./ diff (log(numberOfElements))
pf = polyfit(log(numberOfElements), log(errors), 1);
disp(["EOC:", -pf(1)])
disp("H1")
diffs = diff(-log(errorsH1)) ./ diff (log(numberOfElements))
pf = polyfit(log(numberOfElements), log(errorsH1), 1);
disp(["EOC:", -pf(1)])

subplot(121)
loglog(numberOfElements, errors, '-o')

xlabel('Number of free vertices')
ylabel('||u_{FEM}-u||_{L^2}')
grid on

subplot(122)
loglog(numberOfElements, errorsH1, '-o')

xlabel('Number of free vertices')
ylabel('||u_{FEM}-u||_{H^1}')
grid on
