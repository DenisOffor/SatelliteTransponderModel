#include "processmonitor.h"

#include "processmonitor.h"

#ifdef Q_OS_WIN
#include <windows.h>
#include <psapi.h>
#endif

#include <QtGlobal>

#ifdef Q_OS_WIN
static quint64 fileTimeToUInt64(const FILETIME &ft)
{
    ULARGE_INTEGER li;
    li.LowPart = ft.dwLowDateTime;
    li.HighPart = ft.dwHighDateTime;
    return li.QuadPart;
}
#endif

ProcessMonitor::ProcessMonitor(QObject *parent)
    : QObject(parent)
{
}

void ProcessMonitor::reset()
{
    m_prevProcTime100ns = 0;
    m_prevWallTime100ns = 0;
    m_hasPrevSample = false;

    m_cpuUsagePercent = 0.0;
    m_workingSetMB = 0.0;
    m_privateUsageMB = 0.0;
    m_peakWorkingSetMB = 0.0;
}

void ProcessMonitor::update()
{
#ifdef Q_OS_WIN
    HANDLE hProc = GetCurrentProcess();

    // ---------------- CPU процесса ----------------
    FILETIME ftCreate, ftExit, ftKernel, ftUser;
    FILETIME ftNow;

    if (GetProcessTimes(hProc, &ftCreate, &ftExit, &ftKernel, &ftUser))
    {
        GetSystemTimeAsFileTime(&ftNow);

        const quint64 procKernel = fileTimeToUInt64(ftKernel);
        const quint64 procUser   = fileTimeToUInt64(ftUser);
        const quint64 procTotal  = procKernel + procUser;
        const quint64 wallNow    = fileTimeToUInt64(ftNow);

        if (m_hasPrevSample) {
            const quint64 dProc = procTotal - m_prevProcTime100ns;
            const quint64 dWall = wallNow   - m_prevWallTime100ns;

            if (dWall > 0) {
                SYSTEM_INFO si;
                GetSystemInfo(&si);

                const double cpuCount = qMax<quint32>(1u, si.dwNumberOfProcessors);

                // 100% = полная загрузка приложения по всем ядрам суммарно
                double usage = 100.0 * (double(dProc) / double(dWall)) / cpuCount;

                if (usage < 0.0) usage = 0.0;
                if (usage > 100.0) usage = 100.0;

                m_cpuUsagePercent = usage;
            }
        }

        m_prevProcTime100ns = procTotal;
        m_prevWallTime100ns = wallNow;
        m_hasPrevSample = true;
    }

    // ---------------- RAM процесса ----------------
    PROCESS_MEMORY_COUNTERS_EX pmc;
    if (GetProcessMemoryInfo(hProc,
                             reinterpret_cast<PROCESS_MEMORY_COUNTERS*>(&pmc),
                             sizeof(pmc)))
    {
        constexpr double bytesToMB = 1024.0 * 1024.0;

        m_workingSetMB     = double(pmc.WorkingSetSize) / bytesToMB;
        m_peakWorkingSetMB = double(pmc.PeakWorkingSetSize) / bytesToMB;
        m_privateUsageMB   = double(pmc.PrivateUsage) / bytesToMB;
    }
#else
    m_cpuUsagePercent = 0.0;
    m_workingSetMB = 0.0;
    m_privateUsageMB = 0.0;
    m_peakWorkingSetMB = 0.0;
#endif
}

double ProcessMonitor::cpuUsagePercent() const
{
    return m_cpuUsagePercent;
}

double ProcessMonitor::workingSetMB() const
{
    return m_workingSetMB;
}

double ProcessMonitor::privateUsageMB() const
{
    return m_privateUsageMB;
}

double ProcessMonitor::peakWorkingSetMB() const
{
    return m_peakWorkingSetMB;
}

QString ProcessMonitor::cpuText(int precision) const
{
    return QString("CPU: %1 %").arg(m_cpuUsagePercent, 0, 'f', precision);
}

QString ProcessMonitor::ramText(int precision) const
{
    return QString("RAM: WS %1 MB | Peak %3 MB")
    .arg(m_workingSetMB, 0, 'f', precision)
        .arg(m_peakWorkingSetMB, 0, 'f', precision);
}
