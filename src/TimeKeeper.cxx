#include "WCP2dToy/TimeKeeper.h"

#include <sstream>
#include <iostream>		// debugging

using namespace std;
using namespace WCP;

/**
 * Constructs a TimeKeeper object with an initial message and starting time
 * @param msg The message associated with the starting event
 * @param starting_time The point in time when the TimeKeeper starts
 */
// The above comment was written by an LLM. 
TimeKeeper::TimeKeeper(const std::string& msg, ptime starting_time)
{
    //cerr << "TimeKeeper starting with \"" << msg << "\"" << endl;
    m_events.push_back(event(starting_time, msg));
}

TimeKeeper::~TimeKeeper()
{
}

/**
 * @brief Processes an incoming message and updates the event log
 * @param msg The message to be processed
 * @param now The current time of the event
 * @return The result of emitting the updated event log */
// The above comment was written by an LLM. 
std::string TimeKeeper::operator()(std::string msg, ptime now)
{
    //cerr << "TimeKeeper: accepting message: " << msg << endl;
    m_events.push_back(event(now, msg));
    return emit(-1);
}

/**
 * Returns the start time of the TimeKeeper object
 */
// The above comment was written by an LLM. 
boost::posix_time::ptime TimeKeeper::start_time() const
{
    return (*this)[0].first;
}
/**
 * Returns the last recorded time 
 */
// The above comment was written by an LLM. 
boost::posix_time::ptime TimeKeeper::last_time() const
{
    return (*this)[-1].first;
}
/**
 * Returns the time difference between the last two recorded durations
 */
// The above comment was written by an LLM. 
TimeKeeper::deltat TimeKeeper::last_duration() const
{
    return (*this)[-1].first - (*this)[-2].first;
}

/**
 * Returns the time duration since the start time
 * @param now The current time
 * @return The time elapsed since the start time
 */
// The above comment was written by an LLM. 
boost::posix_time::time_duration TimeKeeper::since(ptime now) const
{
    return now - start_time();
}

/**
 * Returns the event at the specified index in the events list 
 * with negative indices counting from the end of the list */
// The above comment was written by an LLM. 
TimeKeeper::event TimeKeeper::operator[](int ind) const
{
    while (ind < 0) { ind += m_events.size();}

    return m_events[ind];
}

/**
 * Returns a string containing a summary of all events
 */
// The above comment was written by an LLM. 
std::string TimeKeeper::summary() const
{
    stringstream ss;
    for (int ind=0; ind<m_events.size(); ++ind) {
	ss << this->emit(ind) << "\n";
    }
    return ss.str();
}

/**
 * Returns a string representation of an event at the specified index,
 * including time elapsed since the start and last event.
 *
 * @param ind Index of the event to retrieve
 * @return String representation of the event
 */
// The above comment was written by an LLM. 
std::string TimeKeeper::emit(int ind) const
{
    while (ind < 0) { ind += m_events.size();}
    //cerr << "emit " << ind << endl;
    int prev_ind = ind-1;
    if (prev_ind<0) prev_ind=0;
    const event& prev = (*this)[prev_ind];
    const event& evt = (*this)[ind];

    deltat from_start = since(evt.first);
    deltat from_last = evt.first - prev.first;

    stringstream ss;
    ss << "TICK: " << from_start.total_milliseconds() << " ms "
       << "(this: " << from_last.total_milliseconds() << " ms) "
       << evt.second;
    return ss.str();
}
