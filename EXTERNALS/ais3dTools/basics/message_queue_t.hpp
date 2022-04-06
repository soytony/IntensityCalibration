
#include "message_queue_t.h"
#include <math.h>

template <typename ElementType>
MessageQueueT<ElementType>::MessageQueueT()
{
}

template <typename ElementType>
MessageQueueT<ElementType>::~MessageQueueT()
{
}

template <typename ElementType>
ElementType MessageQueueT<ElementType>::findClosestData(double timestamp) const
{
  if (_buffer.rbegin()->first < timestamp)
    return _buffer.rbegin()->second;
  if (_buffer.begin()->first > timestamp)
    return _buffer.begin()->second;

  typename Buffer::const_iterator ub = _buffer.upper_bound(timestamp);
  typename Buffer::const_iterator lb = ub;
  --lb;
  if (fabs(lb->first - timestamp) < fabs(ub->first - timestamp))
    return lb->second;
  else
    return ub->second;
}

template <typename ElementType>
ElementType MessageQueueT<ElementType>::before(double timestamp) const
{
  if (_buffer.size() == 0 || _buffer.begin()->first > timestamp)
    return 0;
  typename Buffer::const_iterator lb = _buffer.upper_bound(timestamp);
  --lb; // now it's the lower bound
  return lb->second;
}

template <typename ElementType>
ElementType MessageQueueT<ElementType>::after(double timestamp) const
{
  if (_buffer.size() == 0 || _buffer.rbegin()->first < timestamp)
    return 0;
  typename Buffer::const_iterator ub = _buffer.upper_bound(timestamp);
  if (ub == _buffer.end())
    return 0;
  return ub->second;
}

template <typename ElementType>
ElementType MessageQueueT<ElementType>::at(double timestamp)
{
  return _buffer[timestamp];
}

template <typename ElementType>
void MessageQueueT<ElementType>::add(double timestamp, ElementType& rd)
{
  _buffer[timestamp] = rd;
}


