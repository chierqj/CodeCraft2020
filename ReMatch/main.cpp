#include <fcntl.h>
#include <sys/mman.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <numeric>
#include <queue>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#define U32 uint32_t

#ifdef LOCAL
#define TRAIN "./data/big/test_data.txt"
#define RESULT "./data/big/result.txt"
#else
#define TRAIN "/data/test_data.txt"
#define RESULT "/projects/student/result.txt"
#endif

const int MAXN = 2000000;
const int NTHREAD = 4;
struct NumBuffer {
  char str[10];
  U32 len;
} __attribute__((packed)) m_MapID[MAXN];
struct ThData {
  U32 answers;
  char Reachable[MAXN];
  std::vector<U32> ReachablePoint;
} ThreadData[NTHREAD];

class XJBG {
 public:
  void Simulation();

 private:
  void ParseInteger(const U32 &x);
  void HandleLoadData(const char *buffer, int pid, U32 st, U32 ed);
  void LoadData();
  void BackSearch(ThData &Data, const U32 &job);
  void ForwardSearch(ThData &Data, const U32 &job);
  void GetNextJob(U32 &job);
  void FindCircle(int pid);
  void SaveAnswer();

 private:
  std::atomic_flag m_lock = ATOMIC_FLAG_INIT;
  U32 m_MaxID = 0;
  U32 m_answers = 0;
  U32 m_jobcur = 0;
  std::vector<U32> m_THIDDom[NTHREAD];
  std::vector<std::tuple<U32, U32, U32>> m_ThEdge[NTHREAD];
  std::vector<U32> m_IDDom;
  std::vector<U32> m_Jobs;
  std::vector<std::pair<U32, U32>> m_Children[MAXN];
  std::vector<U32> m_Parents[MAXN];
  std::vector<std::vector<int>> m_Answer[MAXN];
};

void XJBG::ParseInteger(const U32 &x) {
  U32 num = m_IDDom[x];
  auto &mpid = m_MapID[x];
  if (num == 0) {
    mpid.str[0] = '0';
    mpid.str[1] = ',';
    mpid.len = 2;
  } else {
    char tmp[10];
    int idx = 10;
    tmp[--idx] = ',';
    while (num) {
      tmp[--idx] = num % 10 + '0';
      num /= 10;
    }
    memcpy(mpid.str, tmp + idx, 10 - idx);
    mpid.len = 10 - idx;
  }
}
void XJBG::HandleLoadData(const char *buffer, int pid, U32 st, U32 sz) {
  buffer += st;
  const char *ptr = buffer;
  U32 u = 0, v = 0, w = 0;
  auto &edge = m_ThEdge[pid];
  auto &dom = m_THIDDom[pid];
  while (ptr - buffer < sz) {
    while (*ptr != ',') {
      u = u * 10 + *ptr - '0';
      ++ptr;
    }
    ++ptr;
    while (*ptr != ',') {
      v = v * 10 + *ptr - '0';
      ++ptr;
    }
    ++ptr;
    while (*ptr != '\n') {
      w = w * 10 + *ptr - '0';
      ++ptr;
    }
    ++ptr;
    edge.emplace_back(std::make_tuple(u, v, w));
    dom.emplace_back(u);
    dom.emplace_back(v);
    u = v = w = 0;
  }
  std::sort(dom.begin(), dom.end());
}
void XJBG::LoadData() {
  U32 fd = open(TRAIN, O_RDONLY);
  U32 bufsize = lseek(fd, 0, SEEK_END);
  char *buffer = (char *)mmap(NULL, bufsize, PROT_READ, MAP_PRIVATE, fd, 0);
  close(fd);
  U32 st = 0, block = bufsize / NTHREAD;
  std::thread Th[NTHREAD];
  for (int i = 0; i < NTHREAD; ++i) {
    U32 ed = (i == NTHREAD - 1 ? bufsize : st + block);
    while (buffer[ed] != '\n') ++ed;
    U32 sz = ed - st + 1;
    Th[i] = std::thread(&XJBG::HandleLoadData, this, buffer, i, st, sz);
    st = ed + 1;
  }
  for (int i = 0; i < NTHREAD; ++i) Th[i].join();
  for (int i = 0; i < NTHREAD; ++i) {
    m_IDDom.insert(m_IDDom.end(), m_THIDDom[i].begin(), m_THIDDom[i].end());
  }
  std::sort(m_IDDom.begin(), m_IDDom.end());
  m_MaxID = std::unique(m_IDDom.begin(), m_IDDom.end()) - m_IDDom.begin();
  std::unordered_map<U32, U32> mp;
  for (int i = 0; i < m_MaxID; ++i) {
    mp[m_IDDom[i]] = i;
  }
  std::vector<std::tuple<U32, U32, U32>> E[NTHREAD];
  auto foo = [&](int pid) {
    for (auto &e : m_ThEdge[pid]) {
      const U32 &eu = mp[std::get<0>(e)];
      const U32 &ev = mp[std::get<1>(e)];
      const U32 &ew = std::get<2>(e);
      E[pid].emplace_back(std::make_tuple(eu, ev, ew));
    }
  };
  for (int i = 0; i < NTHREAD; ++i) {
    Th[i] = std::thread(foo, i);
  }
  for (int i = 0; i < NTHREAD; ++i) Th[i].join();
  for (int i = 0; i < NTHREAD; ++i) {
    for (auto &e : E[i]) {
      const U32 &eu = std::get<0>(e);
      const U32 &ev = std::get<1>(e);
      const U32 &ew = std::get<2>(e);
      m_Children[eu].emplace_back(std::make_pair(ev, ew));
      m_Parents[ev].emplace_back(eu);
    }
  }
  for (int i = 0; i < m_MaxID; ++i) {
    if (!m_Children[i].empty() && !m_Children[i].empty()) {
      m_Jobs.emplace_back(i);
      ParseInteger(i);
      std::sort(m_Children[i].begin(), m_Children[i].end());
    }
  }
#ifdef LOCAL
  U32 edges = 0;
  for (int i = 0; i < NTHREAD; ++i) edges += m_ThEdge[i].size();
  std::cerr << "@ u: " << m_MaxID << ", e: " << edges << "\n";
#endif
}
void XJBG::BackSearch(ThData &Data, const U32 &job) {
  for (auto &v : Data.ReachablePoint) Data.Reachable[v] = 0;
  Data.ReachablePoint.clear();
  Data.ReachablePoint.emplace_back(job);
  Data.Reachable[job] = 7;
  for (auto &it1 : m_Parents[job]) {
    if (it1 <= job) continue;
    Data.Reachable[it1] = 7;
    Data.ReachablePoint.emplace_back(it1);
    for (auto &it2 : m_Parents[it1]) {
      if (it2 <= job) continue;
      Data.Reachable[it2] |= 6;
      Data.ReachablePoint.emplace_back(it2);
      for (auto &it3 : m_Parents[it2]) {
        if (it3 <= job || it3 == it1) continue;
        Data.Reachable[it3] |= 4;
        Data.ReachablePoint.emplace_back(it3);
      }
    }
  }
}
void XJBG::ForwardSearch(ThData &Data, const U32 &job) {
  auto &ans = m_Answer[job];
  for (auto &it1 : m_Children[job]) {
    if (it1.first < job) continue;
    for (auto &it2 : m_Children[it1.first]) {
      if (it2.first <= job) continue;
      for (auto &it3 : m_Children[it2.first]) {
        if (it3.first < job || it3.first == it1.first) continue;
        if (it3.first == job) {
          ++Data.answers;
          ans[0].emplace_back(it1.first);
          ans[0].emplace_back(it2.first);
          continue;
        }
        for (auto &it4 : m_Children[it3.first]) {
          if (!(Data.Reachable[it4.first] & 4)) continue;
          if (it4.first == job) {
            ++Data.answers;
            ans[1].emplace_back(it1.first);
            ans[1].emplace_back(it2.first);
            ans[1].emplace_back(it3.first);
            continue;
          }
          if (it4.first == it1.first || it4.first == it2.first) {
            continue;
          }
          for (auto &it5 : m_Children[it4.first]) {
            if (!(Data.Reachable[it5.first] & 2)) continue;
            if (it5.first == job) {
              ++Data.answers;
              ans[2].emplace_back(it1.first);
              ans[2].emplace_back(it2.first);
              ans[2].emplace_back(it3.first);
              ans[2].emplace_back(it4.first);
              continue;
            }
            if (it5.first == it1.first || it5.first == it2.first ||
                it5.first == it3.first) {
              continue;
            }
            for (auto &it6 : m_Children[it5.first]) {
              if (!(Data.Reachable[it6.first] & 1)) continue;
              if (it6.first == job) {
                ++Data.answers;
                ans[3].emplace_back(it1.first);
                ans[3].emplace_back(it2.first);
                ans[3].emplace_back(it3.first);
                ans[3].emplace_back(it4.first);
                ans[3].emplace_back(it5.first);
                continue;
              }
              if (it6.first == it1.first || it6.first == it2.first ||
                  it6.first == it3.first || it6.first == it4.first) {
                continue;
              }
              ++Data.answers;
              ans[4].emplace_back(it1.first);
              ans[4].emplace_back(it2.first);
              ans[4].emplace_back(it3.first);
              ans[4].emplace_back(it4.first);
              ans[4].emplace_back(it5.first);
              ans[4].emplace_back(it6.first);
            }
          }
        }
      }
    }
  }
}

void XJBG::GetNextJob(U32 &job) {
  while (m_lock.test_and_set())
    ;
  if (m_jobcur < m_Jobs.size()) {
    job = m_Jobs[m_jobcur++];
  } else {
    job = -1;
  }
  m_lock.clear();
}

void XJBG::FindCircle(int pid) {
  U32 job = 0;
  auto &Data = ThreadData[pid];
  while (true) {
    GetNextJob(job);
    if (job == -1) break;
    m_Answer[job].resize(5);
    BackSearch(Data, job);
    ForwardSearch(Data, job);
  }
}

void XJBG::SaveAnswer() {
  FILE *fp = fopen(RESULT, "w");
  char tmp[10];
  int tidx = 10;
  U32 tol = m_answers;
  tmp[--tidx] = '\n';
  while (tol) {
    tmp[--tidx] = tol % 10 + '0';
    tol /= 10;
  }
  fwrite(tmp + tidx, 1, 10 - tidx, fp);
  char line[1];
  line[0] = '\n';
  for (int sz = 0; sz < 5; ++sz) {
    for (auto &job : m_Jobs) {
      if (m_Answer[job].empty()) continue;
      const auto &ans = m_Answer[job][sz];
      if (ans.empty()) continue;
      const auto &mpjob = m_MapID[job];
      int idx = 0;
      for (auto &v : ans) {
        if (!idx) fwrite(mpjob.str, 1, mpjob.len, fp);
        ++idx;
        if (idx == sz + 2) {
          idx = 0;
          fwrite(m_MapID[v].str, 1, m_MapID[v].len - 1, fp);
          fwrite(line, 1, 1, fp);
        } else {
          fwrite(m_MapID[v].str, 1, m_MapID[v].len, fp);
        }
      }
    }
  }
  fclose(fp);
}
void XJBG::Simulation() {
  LoadData();
  /*
  std::thread Th[NTHREAD];
  for (int i = 0; i < NTHREAD; ++i) {
    Th[i] = std::thread(&XJBG::FindCircle, this, i);
  }
  for (auto &it : Th) it.join();
  for (auto &it : ThreadData) m_answers += it.answers;
  SaveAnswer();
  */

#ifdef LOCAL
  std::cerr << "@ answers: " << m_answers << "\n";
  for (int i = 0; i < NTHREAD; ++i) {
    std::cerr << "@ " << i << ": " << ThreadData[i].answers << "\n";
  }
#endif
}
int main() {
  XJBG *xjbg = new XJBG();
  xjbg->Simulation();
  return 0;
}