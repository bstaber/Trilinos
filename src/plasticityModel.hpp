/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef PLASTICITYMODEL_HPP
#define PLASTICITYMODEL_HPP

class subdifferentialTresca
{
  subdifferentialTresca();
  ~subdifferentialTresca();

  void yield_function();
  void integrate();
  void newton_left();
  void newton_right();
  void newton_smooth();
};
#endif
