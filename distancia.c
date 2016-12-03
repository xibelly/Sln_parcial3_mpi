//funcion que calcula la distancia entre las particulas
double distance(double xii, double yii, double zii, double xjj, double yjj, double zjj)
{

  double dist, dx, dy, dz;
  
  dx = (xii-xjj)*(xii-xjj);
  dy = (yii-yjj)*(yii-yjj);
  dz = (zii-zjj)*(zii-zjj);
  
  dist = sqrt(dx + dy + dz); 
  
  return dist;

}
