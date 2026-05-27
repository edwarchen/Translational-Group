# Translational Group

转化组仓库。**main 分支为团队入口，项目代码在个人分支上管理。**

## 分支模型

```
main  (只读 — 仅含本 README + CI 配置)
│
├── edwarchen/dev  ───  projects/TCR_Pipeline_Optim/      # 成员项目
├── zhangsan/dev   ───  projects/his-project/             # 成员项目
├── lisi/dev       ───  projects/her-project/             # 成员项目
└── ...
```

- **main 分支不含项目代码**，仅保留团队入口和 CI 配置
- 每个成员在自己的 `{username}/{suffix}` 分支上独立开发
- `projects/` 目录下各自管理项目，互不干扰
- 合并到 main 通过 Pull Request + 1 人审批

## 当前成员与项目

| 成员 | 分支 | 项目 |
|---|---|---|
| [@edwarchen](https://github.com/edwarchen) | [`edwarchen/dev`](https://github.com/edwarchen/Translational-Group/tree/edwarchen/dev) | TCR 管线优化（MiXCR 预设对比、VDJTools 交叉验证） |

## 新成员加入

### 1. 获得仓库访问权限
管理员在 **Settings → Collaborators → Add people** 邀请你的 GitHub 账号。

### 2. 创建你的分支
接受邀请后，在仓库 **Actions → Create User Branch → Run workflow** 输入你的 GitHub 用户名，自动创建 `{username}/dev` 分支。

### 3. 克隆并开始工作

```bash
git clone https://github.com/edwarchen/Translational-Group.git
cd Translational-Group
git fetch origin
git checkout your-username/dev

# 创建你的项目目录
mkdir -p projects/your-project
cd projects/your-project
# ... 开始工作 ...
```

### 4. 提交更改

```bash
git add projects/your-project/
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
