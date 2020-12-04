

function splineInter(a, b, n, f) {

  let h;



}

function getH(a, b, n) {
  return (b-a)/n;
}

function getX(a, b, n, h) {
  let x = [];

  for (let i=0;i<n;i++) {
    x[i] = a + i*h;
  }
  x[n] = b;

  return x;
}

function getS3() {

}

//Способ 1 (упрощенный)
function getM1(f, n, h) {
  let m;
  m[0] = (4*f[1] - f[2] - 3*f[0])/(2*h);

  for (let i=0; i<n; i++) {
    m[i] = (f[i+1] - f[i-1])/(2*h);
  }

  m[n] = (3*f[n] - f[n-2] - 3*f[n-1])/(2*h);
  return m;
}

//Способ 2
function getM2(f2, n, h) {

  for (let i = 0;i<n+1;i++) {
    m[i] = f2[i];
  }

  return m;
}

//Способ 3 (глобальный)

function getS3Plus3() {

}

function getS3Minus3() {

}

function mmm() {

}

function* getMFirstLast(type, f, n, h, m1, mnMinus1, f20, f2n) {
  while (true) {
    let m0, mn;

    if (type = 0) {
      m0 = f[0];
      mn = f[n];
    }
    if (type = 0) {
      m0 = (1/6*h)*(-11*f[0] + 18*f[1] - 9*f[2] + 2*f[3]);
      mn = (1/6*h)*(11*f[n] - 18*f[n-1] + 9*f[n-2] - 2*f[n-3]);
    }
    if (type = 0) {
      m0 = (-m1/2) + ( (3/2)*((*f[1] - f[0])/h) ) + (h/4*f20)
      mn = (-mnMinus1/2) + ( (3/2)*((*f[n] - f[n-1])/h) ) + (h/4*f2n)
    }

    yield m0;
    yield mn;
  }
}
