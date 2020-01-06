#ifndef ERROR_MESSAGER_HPP
#define ERROR_MESSAGER_HPP

#include <exception>

class NotImplementedException : public std::logic_error
{
public:
  NotImplementedException() : std::logic_error("Function not yet implemented") { };
};

#endif
