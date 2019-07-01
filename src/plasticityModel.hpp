/*
Brian Staber (brian.staber@gmail.com)
*/

#ifndef PLASTICITYMODEL_HPP
#define PLASTICITYMODEL_HPP

class plasticityModel
{
  plasticityModel();
  ~plasticityModel();

  virtual void yield_function() = 0;
  virtual void hardening_function();
  virtual void return_mapping() = 0;
  void plasticMultiplier();
};
#endif
