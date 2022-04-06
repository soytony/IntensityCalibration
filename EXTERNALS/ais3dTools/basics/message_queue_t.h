#ifndef MESSAGE_QUEUE_T_H
#define MESSAGE_QUEUE_T_H

#include <map>

namespace Ais3dTools{
  
template <typename ElementType>
class MessageQueueT
{
  public:
    typedef std::map<double, ElementType>           Buffer;

  public:
    MessageQueueT();
    ~MessageQueueT();

    void add(double timestamp, ElementType& rd);

    ElementType findClosestData(double timestamp) const;

    ElementType before(double timestamp) const;
    ElementType after(double timestamp) const;

    ElementType at(double timestamp);
	
    const Buffer& buffer() const {return _buffer;}

  protected:
    Buffer _buffer;
};

#include "message_queue_t.hpp"

} //end namespace

#endif

