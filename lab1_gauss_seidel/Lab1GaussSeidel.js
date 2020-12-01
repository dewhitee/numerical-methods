function gaussSeidel (a,b,e) {
  let xNew = new Array();
  let x = new Array();
  let stop = false;
  let itter = 0;
  let eRes = [];
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
          xNew[i] += -a[i][j] * x[j];
        }
        if ((i != j) && (j<i)) {
          xNew[i] += -a[i][j] * xNew[j];
        }
      }
      xNew[i] += b[i];
      xNew[i] *= 1/a[i][i];

    }

    //проверка e
    stop = false;
    for (let i = 0; i<x.length;i++) {

      if (Math.abs(xNew[i] - x[i]) <= e) {
        stop = true;
      } else {
        stop = false;
        break;
      }
    }
    //выход из while
    if (stop) {
      for (let i = 0; i<x.length;i++) {
        eRes.push((xNew[i] - x[i]).toFixed(6));
        xNew[i] = xNew[i].toFixed(6);
      }
      break;
    }

    //выход если достигли лимит итераций
    if (itter>99) {
      for (let i = 0; i<x.length;i++) {
        eRes.push((xNew[i] - x[i]).toFixed(6));
        xNew[i] = xNew[i].toFixed(6);
      }
      return "\n---------\n\nДостигнут лимит итераций\ne = " + e + "\nXi+1 - Xi = ["+ eRes +"]\nX =  [" + xNew + "]\n\n---------";
    }

    //xNew становятьося x
    for (let i = 0; i<x.length; i++) {
      x[i] = xNew[i];
    }
    itter++;

  }
  //console.log("X = [" + xNew + "]");

  return "\n---------\n\ne = " + e + "\nXi+1 - Xi = ["+ eRes +"]\nX =  [" + xNew + "]\n\n---------";
}
