# no_ecoqos.ps1 -- keeps Windows from running the simulation's R processes in
# efficiency mode (EcoQoS), which on this Core Ultra laptop parks background work on
# the slow cores (2026-10-17: workers at ~11% of a core, replicates 3x slower).
# Clears the per-process execution-speed throttling flag on every Rscript/Rterm
# process; with -Watch it repeats every 15 s, since run_scenario starts fresh
# workers for every chunk. Changes no system setting.
param([switch]$Watch)
Add-Type @"
using System; using System.Runtime.InteropServices;
public static class Eco {
  [StructLayout(LayoutKind.Sequential)] public struct PPTS { public uint Version, ControlMask, StateMask; }
  [DllImport("kernel32.dll", SetLastError=true)] public static extern IntPtr OpenProcess(uint a, bool i, int pid);
  [DllImport("kernel32.dll", SetLastError=true)] public static extern bool SetProcessInformation(IntPtr h, int cls, ref PPTS s, int size);
  [DllImport("kernel32.dll")] public static extern bool CloseHandle(IntPtr h);
  public static bool Off(int pid) {
    IntPtr h = OpenProcess(0x0200, false, pid);            // PROCESS_SET_INFORMATION
    if (h == IntPtr.Zero) return false;
    var s = new PPTS { Version = 1, ControlMask = 1, StateMask = 0 };  // execution speed: not throttled
    bool ok = SetProcessInformation(h, 4, ref s, Marshal.SizeOf(s));   // ProcessPowerThrottling
    CloseHandle(h); return ok;
  }
}
"@
$done = @{}
do {
  foreach ($p in Get-Process Rscript, Rterm -ErrorAction SilentlyContinue) {
    if (-not $done.ContainsKey($p.Id)) { $done[$p.Id] = [Eco]::Off($p.Id) }
  }
  if (-not $Watch) { $done.GetEnumerator() | ForEach-Object { "{0} {1}" -f $_.Key, $_.Value } }
  if ($Watch) { Start-Sleep 15 }
} while ($Watch -and (Get-Process Rscript, Rterm -ErrorAction SilentlyContinue))
