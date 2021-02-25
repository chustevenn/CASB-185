plot(log(out.simout.time), out.simout.signals.values, 'LineWidth',1)
xlim([0,7])
xlabel("ln(t)")
ylim([0,210])
ylabel("N(t)")
legend(["N_r(t) - Recovered","N_s(t) - Susceptible","N_i(t) - Infected"],"Location","northwest")
grid on
