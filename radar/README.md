# coding-rna 研究雷达（看板）

- 页面：https://claude.ai/artifact/2XzrLqNP88edE6MnnGt3d3（私有；数据库只有 Claude/所有者可写；页面实时订阅）
- 数据模型与 TE 看板相同（collections：projects/crna、nodes、claims、decisions、queues、log、sota），字段说明见 /Users/jiwang/Desktop/TE/TE-GLM-study/00_overview/radar/README.md
- `seed/` 为 2026-10-22 子代理盘点的初始数据（SUMMARY.md 为现状总结）；`crna-radar.html` 为页面源码。
- 维护规则：实验状态变化 → update node；需要用户拍板 → 新建 decision（why/options/consequence/recommended）；用户在页内点选后 decision 会带 acknowledged:false，每次开工先读取并确认；当天的事写进当天 log；更新 projects/crna 的 focus/nextCheckpoint/updatedAt。
- 工作方式：先论文架构（claims × figures），实验只为补证据；新想法先进 decisions 挂到某条论点。
