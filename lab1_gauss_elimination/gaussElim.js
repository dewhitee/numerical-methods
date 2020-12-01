function gaussMain(a,b) {
  let xResult = [];
  let stop = false;
  let positionMain;
  let mainElement;
  //console.log("A = " + a);
  //console.log("B = " + b);

  for (let m = 0;m<a.length; m++) {
    positionMain = mainElemntPosition(a,m);
    mainElement = a[positionMain[0]][positionMain[1]];
    a = swapRowA(a,m,positionMain);
    b = swapRowB(b,m,positionMain);
    b = divBOnMain(b,a,mainElement,m);
    a = divAOnMain(a,mainElement,m);
    //console.log("A = " + a);
    //console.log("B = " + b);
  }

  xResult = findX(a,b);

  return "\n---------\n\nX =  [" + xResult + "]\n\n---------";
}

function mainElemntPosition(a,j) {
  let positionBigI = 0;

  for (let i = j;i<a.length; i++) {
    if (a[positionBigI][j] < a[i][j]) {
      positionBigI = i;
    }
  }

  positionBigEl = [positionBigI,j];

  return positionBigEl;
}

function swapRowA(a,i,positionBigEl) {
  let tmpA = 0;

  for (let j = 0;j<a[i].length; j++) {
    tmpA = a[i][j];
    a[i][j] = a[positionBigEl[0]][j];
    a[positionBigEl[0]][j] = tmpA;
  }

  return a;
}

function swapRowB(b,i,positionBigEl) {
  let tmpB = 0;

  tmpB = b[i];
  b[i] = b[positionBigEl[0]];
  b[positionBigEl[0]] = tmpB;

  return b;
}


function divAOnMain(a,mainEl,i) {

  for (let j = 0;j<a[i].length; j++) {
    a[i][j] /= mainEl;
  }

  if (!((i+1)<a.length)) {
    return a;
  }

  for (let k = i+1;k<a.length; k++) {
    for (let j = i+1;j<a[i].length; j++) {
      a[k][j] -= a[i][j]*a[k][i];
    }
  }

  for (let k = i+1;k<a[i].length; k++) {
    a[k][i] = 0;
  }

  return a;
}

function divBOnMain(b,a,mainEl,i) {
    b[i] /= mainEl;

    if (!((i+1)<b.length)) {
      return b;
    }

    for (let m = i+1; m<b.length; m++) {
      b[m] -= b[i]*a[m][i];;
    }

  return b;
}

function findX(a,b) {
  let x = [];

  for (let i = a.length-1;i>-1;i--) {

    x[i] = b[i];
    if (i+1 < a[i].length) {
      for (let j=i+1;j<a[i].length;j++){
        x[i] += -a[i][j]*x[j];
      }
    }
    x[i] /= a[i][i];

  }
  return x;
}
