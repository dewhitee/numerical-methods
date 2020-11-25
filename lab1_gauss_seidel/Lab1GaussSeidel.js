function gaussSeidel (a,b,e) {
  let xNew = new Array();
  let x = new Array();
  let stop = false;
  let itter = 0;
  //Ввод данных
  //x = [x1,x2 ...]
  //b = [b1,b2 ...]
  //a = [[a11,a12 ...],
  //     [a21,a22 ...],
  //     [a31,a32 ...],
  //     [a41,a42 ...]]
  //       ......

  for (let i = 0; i<b.length; i++) {
    x[i] = b[i];
  }

  while (true) {

    //итерация
    for (let i = 0; i<x.length; i++) {

      xNew[i] = 0;
      for (let j = 0; j<a[i].length; j++) {
        if ((i != j) && (j>i)) {
          xNew[i] += a[i][j] * x[j];
        }
        if ((i != j) && (j<i)) {
          xNew[i] += a[i][j] * xNew[j];
        }
        //console.log(xNew[i]);
      }
      xNew[i] += b[i];
      //console.log("x = " + xNew[i]);

    }
    console.log("xNew = " + xNew);

    //проверка e
    stop = false;
    for (let i = 0; i<x.lenght;i++) {

      console.log(Math.abs(xNew[i] - x[i]));
      if (Math.abs(xNew[i] - x[i]) < e) {
        stop = true;
      } else {
        stop = false;
        break;
      }
    }
    //выход из while
    if (stop || (itter>99)) {
      break;
    }

    //xNew становятьося x
    for (let i = 0; i<x.length; i++) {
      x[i] = xNew[i];
    }
    itter++;

  }
  console.log("xNew = " + xNew);

  return "x = { " + xNew + " }";
}
