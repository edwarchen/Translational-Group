# Translational Group

肿瘤转化研究协作仓库。每个成员在自己的分支下管理独立项目。

## 仓库结构

```
Translational-Group/
├── projects/
│   ├── TCR_Pipeline_Optim/     # edwarchen: TCR 分析管线优化
│   └── (your project)/         # 你的项目
└── README.md
```

## 当前项目

| 项目 | 负责人 | 分支 | 说明 |
|---|---|---|---|
| [TCR_Pipeline_Optim](projects/TCR_Pipeline_Optim/) | @edwarchen | `edwarchen/dev` | MiXCR 预设对比与标准化 TCR 分析管线 |

## 协作方式

### 分支策略

```
main (受保护，PR + 审批强制)
├── edwarchen/dev           # 个人分支
├── username/dev            # 你的分支
└── ...
```

- **main** 受保护，所有改动必须通过 Pull Request + 审批
- 每个人在自己的 `{username}/{suffix}` 分支上独立开发
- 互不干扰，各自管理各自的项目

### 新成员加入

1. 管理员在 **Settings → Collaborators** 邀请你
2. 接受邀请后，在 **Actions → Create User Branch** 输入你的 GitHub 用户名，自动创建 `{username}/dev` 分支
3. 克隆并切换到自己的分支：

```bash
git clone https://github.com/edwarchen/Translational-Group.git
cd Translational-Group
git fetch origin
git checkout your-username/dev
```

4. 在 `projects/` 下创建自己的项目目录，开始工作

### 提交与合并

```bash
git add projects/your-project/
git commit -m "feat: description"
git push

# 在 GitHub 网页创建 Pull Request: your-username/dev → main
# 需要至少 1 人审批
```

## 分支保护规则

- **Require a pull request before merging** ✓
- **Require 1 approval** ✓
- **Allow force pushes** ✗
- **Allow deletions** ✗
