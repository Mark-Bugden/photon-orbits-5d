"""Check the grouped closed form of the quartic F(x) = P'Delta - P Delta' quoted in the report,
and the boundary (constant-theta) formulas for Phi^2, Psi^2."""
import sympy as sp
x, mu, a, b, Phi, Psi, K, y0 = sp.symbols('x mu a b Phi Psi K y0', real=True)
P = x**3 + (a**2 + b**2)*x**2 + ((a**2 - b**2)*(Phi**2 - Psi**2) + a**2*b**2 + mu*(a**2 + b**2 + 2*a*Phi + 2*b*Psi))*x \
    + mu*(a*b + b*Phi + a*Psi)**2 + (a**2 - b**2)*(b**2*Phi**2 - a**2*Psi**2)
D = x**2 + (a**2 + b**2 - mu)*x + a**2*b**2
F = sp.expand(sp.diff(P, x)*D - P*sp.diff(D, x))
c3 = 2*(a**2 + b**2 - mu)
c2 = (b**2 - a**2)*(Phi**2 - Psi**2) - 2*mu*(a*Phi + b*Psi) + a**4 + 4*a**2*b**2 + b**4 - 2*mu*(a**2 + b**2)
c1 = 2*(b**2*(b**2 - a**2 - mu)*Phi**2 + a**2*(a**2 - b**2 - mu)*Psi**2 - 2*mu*a*b*(Phi*Psi + b*Phi + a*Psi) + a**2*b**2*(a**2 + b**2 - mu))
c0 = b**2*((b**2 - mu)**2 - a**2*b**2)*Phi**2 + a**2*((a**2 - mu)**2 - a**2*b**2)*Psi**2 + 2*a*b*mu*(mu - a**2 - b**2)*Phi*Psi \
     + 2*a*b**2*mu*(mu - b**2)*Phi + 2*a**2*b*mu*(mu - a**2)*Psi + a**2*b**2*(a**2*b**2 + mu**2)
print("grouped F matches:", sp.expand(F - (x**4 + c3*x**3 + c2*x**2 + c1*x + c0)) == 0)
# chi = P - K*Delta reproduces calX (E=1) -- already checked in check_symbolic.py; here check boundary formulas:
Theta = lambda yy: a**2*(1 - yy) + b**2*yy + K - Phi**2/yy - Psi**2/(1 - yy)
Phi2 = y0**2*(K + 2*a**2*(1 - y0) + b**2*(2*y0 - 1))
Psi2 = (1 - y0)**2*(K + 2*b**2*y0 + a**2*(1 - 2*y0))
yy = sp.symbols('yy')
T = Theta(yy).subs({Phi**2: Phi2, Psi**2: Psi2})
T = a**2*(1 - yy) + b**2*yy + K - Phi2/yy - Psi2/(1 - yy)
print("Theta(y0) = 0 :", sp.simplify(T.subs(yy, y0)) == 0)
print("Theta'(y0) = 0:", sp.simplify(sp.diff(T, yy).subs(yy, y0)) == 0)
