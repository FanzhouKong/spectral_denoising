import numpy as np


class Run(dict):
    def __init__(self, run=None):
        super(Run, self).__init__()
        self["rt_list_ms1"] = []

        self["scan_list_ms1"] = []
        self["scan_list_ms2"] = []

        if run is not None:
            self.update(run)

    # Get scans
    def iter_ms1(self):
        for ms1_scan in self["scan_list_ms1"]:
            yield ms1_scan

    def iter_ms2(self):
        for ms2_scan in self["scan_list_ms2"]:
            yield ms2_scan

    def get_nearest_ms1_scan_index(self, rt, direction="left"):
        if direction == "left":
            scan_idx = np.searchsorted(self["rt_list_ms1"], rt, side="right") - 1
        elif direction == "right":
            scan_idx = np.searchsorted(self["rt_list_ms1"], rt, side="left") + 1
        if scan_idx < 0:
            scan_idx = 0
        elif scan_idx >= self.get_ms1_scan_num():
            scan_idx = self.get_ms1_scan_num() - 1
        return scan_idx

    def get_ms1_scan(self, index):
        return self["scan_list_ms1"][index]

    def get_ms2_scan(self, index):
        return self["scan_list_ms2"][index]

    def get_ms1_scan_num(self):
        return len(self["scan_list_ms1"])

    def get_ms2_scan_num(self):
        return len(self["scan_list_ms2"])

    # Adding data
    def add_scan(self, scan):
        scan = Scan(scan)
        if scan["_ms_level"] == 1:
            self["scan_list_ms1"].append(scan)
            self["rt_list_ms1"].append(scan["rt"])
        elif scan["_ms_level"] == 2:
            self["scan_list_ms2"].append(scan)

    def finish_adding_data(self):
        self["rt_list_ms1"] = np.array(self["rt_list_ms1"])


class Scan(dict):
    def __init__(self, scan=None):
        super(Scan, self).__init__()
        if scan is not None:
            self.update(scan)

        if "peaks" in self:
            peaks = np.asarray(self["peaks"], dtype=np.float32)
            self["peaks"] = peaks[np.argsort(peaks[:, 0])]

    def get_nearest_mz_index(self, mz, direction="left"):
        if direction == "left":
            mz_idx = np.searchsorted(self["peaks"][:, 0], mz, side="right") - 1
        elif direction == "right":
            mz_idx = np.searchsorted(self["peaks"][:, 0], mz, side="left")
        if mz_idx < 0:
            mz_idx = 0
        elif mz_idx >= len(self["peaks"]):
            mz_idx = len(self["peaks"]) - 1
        return mz_idx

    def get_intensity_between_mz(self, mz_left, mz_right):
        mz_idx_left = self.get_nearest_mz_index(mz_left, direction="right")
        mz_idx_right = self.get_nearest_mz_index(mz_right, direction="left") + 1
        intensity = np.sum(self["peaks"][mz_idx_left:mz_idx_right, 1])
        return intensity
