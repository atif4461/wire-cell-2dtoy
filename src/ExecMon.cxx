#include "WCP2dToy/ExecMon.h"
#include <sstream>

using namespace WCP;

ExecMon::ExecMon(const std::string& msg, TimeKeeper::ptime starting_time)
    : tk(msg, starting_time)
    , mu(msg) { }
    
ExecMon::~ExecMon() { }

/**
 * @brief Executes monitoring operation with message time and memory usage
 * @param msg input message
 * @param now current time
 * @param mumu memory usage data
 * @return formatted string containing time and memory information
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
 * Returns a string containing time and memory usage summaries.
 */
// The above comment was written by an LLM. 
std::string ExecMon::summary() const
{
    std::stringstream ss;
    ss << "Time summary:\n" << tk.summary() << "\nMemory usage:\n" << mu.summary();
    return ss.str();
}

