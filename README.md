# 转化组

转化组协作仓库。个人在自己的分支上开发，项目成熟并经审查后合并入 main。

## 分支结构

```
main  (受保护 — 团队入口 + 已审查通过的成熟项目)
│
├── edwarchen/dev  ───  TCR_Pipeline_Optim/              # 开发中
├── zhangsan/dev   ───  his-project/                     # 开发中
├── lisi/dev       ───  her-project/                     # 开发中
└── ...
```

- **main 分支**：存放团队入口文档、CI 配置，以及**经 PR 审查通过、可融入主 pipeline 的成熟项目代码**
- **个人分支**：每个成员在 `{username}/{suffix}` 上自由开发，不受限制
- **合并流程**：个人分支 → PR → 审查 → merge 到 main

## 当前项目

| 项目 | 负责人 | 状态 | 分支 |
|---|---|---|---|
| TCR 管线优化 | [@edwarchen](https://github.com/edwarchen) | 开发中 | [`edwarchen/dev`](https://github.com/edwarchen/Translational-Group/tree/edwarchen/dev) |

*main 暂无已合并的成熟项目。*

## 新成员加入

### 1. 获得仓库访问权限
管理员在 **Settings → Collaborators → Add people** 邀请你的 GitHub 账号。

### 2. 创建你的分支
接受邀请后，在仓库 **Actions → Create User Branch → Run workflow** 输入你的 GitHub 用户名，自动创建 `{username}/dev` 分支。

### 3. 克隆并开始工作

```bash
git clone https://github.com/edwarchen/Translational-Group.git
cd Translational-Group
git fetch origin  # 获取远程仓库最新的状态
git checkout your-username/dev  # 切换到个人分支

# 创建你的项目目录
mkdir your-project
cd your-project
# ... 开始工作 ...
```

### 4. 提交更改

```bash
git add your-project/
git commit -m "feat: description of changes"
git push

# 在 GitHub 网页创建 Pull Request → your-username/dev → main
# 需要至少 1 人审批通过后合并
```

## 分支保护规则

| 规则 | 状态 |
|---|---|
| Require pull request before merging | ✅ |
| Require 1 approval | ✅ |
| Allow force pushes | ❌ |
| Allow deletions | ❌ |
