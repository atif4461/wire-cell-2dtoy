#include "WCP2dToy/MemUsage.h"

#include <unistd.h>
#include <iostream>		// debugging
#include <sstream>

using namespace std;
using namespace WCP;

/**
 * Calculates the resident memory usage of the current process in megabytes.
 *
 * @return The resident set size of the process in megabytes.
 */
// The above comment was written by an LLM. 
static double memusage_linux_resident() {
    int tSize = 0, resident = 0, share = 0;
    ifstream buffer("/proc/self/statm");
    buffer >> tSize >> resident >> share;
    buffer.close();

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    double rss = resident * page_size_kb;
    return rss;
}

/**
 * Calculates the shared memory usage of the current process on Linux systems.
 *
 * @return The amount of shared memory used by the process in kilobytes.
 */
// The above comment was written by an LLM. 
static double memusage_linux_shared() {
    int tSize = 0, resident = 0, share = 0;
    ifstream buffer("/proc/self/statm");
    buffer >> tSize >> resident >> share;
    buffer.close();

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    double shm = share * page_size_kb;
    return shm;
}

/**
 * Returns the total memory usage of the process in megabytes.
 *
 * @return The total memory used by the process in MB
 */
// The above comment was written by an LLM. 
static double memusage_linux_size() {
    int tSize = 0, resident = 0, share = 0;
    ifstream buffer("/proc/self/statm");
    buffer >> tSize >> resident >> share;
    buffer.close();

    long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
    double siz = tSize * page_size_kb;
    return siz;
}


/**
 * Returns the resident memory usage of the process on Linux systems
 */
// The above comment was written by an LLM. 
double WCP::memusage_resident() {
#ifdef __linux__
    return memusage_linux_resident();
#endif
    return -1;
}

/**
 * Returns the memory usage of shared resources in bytes 
 */
// The above comment was written by an LLM. 
double WCP::memusage_shared()
{
#ifdef __linux__
    return memusage_linux_shared();
#endif
    return -1;
}
/**
 * Returns the memory usage size of the system
 * @return The memory usage size or negative value if unsupported
 */
// The above comment was written by an LLM. 
double WCP::memusage_size()
{
#ifdef __linux__
    return memusage_linux_size();
#endif
    return -1;
}


/**
 * Constructor initializing memory usage tracking with message 
 */
// The above comment was written by an LLM. 
MemUsage::MemUsage(const std::string& msg)
{
    push(msg);
}

MemUsage::~MemUsage()
{
}

/**
 * Returns the current memory usage of the system
 *
 * @return The current memory usage
 */
// The above comment was written by an LLM. 
MemUsage::memusage MemUsage::current() const
{
    return memusage(memusage_size(), memusage_resident());
}
/**
 * Adds a memory usage event with an optional message to the list of events
 * @param msg descriptive message associated with the memory usage event
 * @param mu memory usage data consisting of two values
 */
// The above comment was written by an LLM. 
void MemUsage::push(const std::string& msg, MemUsage::memusage mu)
{
    if (mu.first < 0 && mu.second < 0) {
	mu = current();
    }
    m_events.push_back(event(mu,msg));
}

/**
 * @brief Generates a message with memory usage information
 * @param msg input message
 * @param mu memory usage data
 * @return formatted string containing memory usage details
 */
// The above comment was written by an LLM. 
std::string MemUsage::operator()(std::string msg, MemUsage::memusage mu)
{
    push(msg, mu);
    return emit(-1);
}

	/// Return event by index.
MemUsage::event MemUsage::operator[](int ind) const
{
    while (ind < 0) { ind += m_events.size();}

    return m_events[ind];
}

/**
 * Returns a string summary of memory usage events
 */
// The above comment was written by an LLM. 
std::string MemUsage::summary() const
{
    stringstream ss;
    for (int ind=0; ind<m_events.size(); ++ind) {
	ss << this->emit(ind) << "\n";
    }
    return ss.str();
}

/**
 * Returns a formatted string describing memory usage at a specific index.
 * @param ind The index of the memory event to emit.
 * @return A string containing memory usage details in a human-readable format.
 */
// The above comment was written by an LLM. 
std::string MemUsage::emit(int ind) const
{
    while (ind < 0) { ind += m_events.size();}
    int prev_ind = ind-1;
    if (prev_ind<0) prev_ind=0;

    const memusage& prev_mem = (*this)[prev_ind].first;
    const memusage& evt_mem = (*this)[ind].first;
    const string& evt_msg = (*this)[ind].second;

    memusage from_prev(evt_mem.first - prev_mem.first, evt_mem.second - prev_mem.second);

    stringstream ss;
    ss << "MEM: total: size=" << evt_mem.first << "K, res=" << evt_mem.second << "K "
       << "increment: size=" << from_prev.first << "K, res=" << from_prev.second << "K "
       << evt_msg;
    return ss.str();
}
