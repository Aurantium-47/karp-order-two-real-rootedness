"""Exact low-degree identities for Karp's order-two transform.
Requires SymPy and assertions enabled.
"""
import sympy as s

def verify():
    x = s.Symbol("x")
    f=s.symbols('f0:8')
    expected={
        3:6*(f[1]*f[2]-f[0]*f[3])*(x+1),
        4:12*(x+1)*((f[2]**2-f[0]*f[4])*x+2*(f[1]*f[3]-f[0]*f[4])),
        5:20*(x+1)*(x+2)*((2*f[2]*f[3]-f[0]*f[5]-f[1]*f[4])*x+3*(f[1]*f[4]-f[0]*f[5])),
        6:30*(x+1)*(x+2)*((-f[0]*f[6]-2*f[1]*f[5]+f[2]*f[4]+2*f[3]**2)*x**2
            +(-7*f[0]*f[6]-2*f[1]*f[5]+7*f[2]*f[4]+2*f[3]**2)*x+12*(f[1]*f[5]-f[0]*f[6])),
        7:42*(x+1)*(x+2)*(x+3)*((-f[0]*f[7]-3*f[1]*f[6]-f[2]*f[5]+5*f[3]*f[4])*x**2
            +(-9*f[0]*f[7]-7*f[1]*f[6]+11*f[2]*f[5]+5*f[3]*f[4])*x+20*(f[1]*f[6]-f[0]*f[7]))
    }
    for n,want in expected.items():
        actual=sum(s.binomial(n,k)*f[k]*f[n-k]*(s.rf(x,k)*s.rf(x,n-k)-s.rf(x+1,k)*s.rf(x-1,n-k)) for k in range(n+1))
        assert s.expand(actual-want)==0
    p,q,r=s.symbols('p q r')
    delta6=(7*p+2*q)**2-48*r*(p+2*q)
    bound6=s.Rational(1033,25)*p**2+s.Rational(316,25)*p*q+4*q**2
    assert s.expand(delta6.subs(r,4*p/25)-bound6)==0
    L=p+4*q+5*r; M=9*p+16*q+5*r
    assert s.expand(M*M-80*p*L-(L*L+(24*q-64*p)*L+(8*p+12*q)**2))==0
    print("Karp n=3,...,7 symbolic identities verified.")

if __name__ == "__main__":
    verify()
