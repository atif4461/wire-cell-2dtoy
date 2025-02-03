#include "WCP2dToy/ExecMon.h"
#include <sstream>

using namespace WCP;

ExecMon::ExecMon(const std::string& msg, TimeKeeper::ptime starting_time)
    : tk(msg, starting_time)
    , mu(msg) { }
    
ExecMon::~ExecMon() { }

/**
 * @brief Generates a formatted string containing time and memory usage information
 * @param msg input message
 * @param now current time
 * @param mumu memory usage data
 * @return formatted string with time and memory details
 */
// The above comment was written by an LLM. 
std::string ExecMon::operator()(std::string msg, TimeKeeper::ptime now, MemUsage::memusage mumu)
{
    std::stringstream ss;
    ss << "Time: " << tk(msg,now) << "\n"
       << "Memory: " << mu(msg,mumu);
    return ss.str();
}


/**
 * Returns a string containing time and memory usage summaries
 */
// The above comment was written by an LLM. 
std::string ExecMon::summary() const
{
    std::stringstream ss;
    ss << "Time summary:\n" << tk.summary() << "\nMemory usage:\n" << mu.summary();
    return ss.str();
}

