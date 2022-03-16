fplot(@(x) 3.11.*x)
hold on
fplot(@(x)  1.61E-03 + 4.45.*x)
hold off

xlim([0 2e-3])
ylim([0 4e-3])
xlabel('Mole fraction of KOH in water, x_i, (unitless)')
ylabel('Mole fraction of KOH in hydrogen, y_i, (unitless)')
legend('Equilibrium line', 'Operating line')