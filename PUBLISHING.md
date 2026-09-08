# Publish this snapshot to GitHub

The directory is already an initialized local Git repository. It has no remote
and no commit because the GitHub owner, repository URL, and commit identity are
user-specific.

## 1. Create the empty GitHub repository

On GitHub, select **New repository** and use a name such as `flightcollapse`.
For an initial collaborator evaluation, **Private** is the safer default. Do not
pre-populate it with a README, `.gitignore`, or license because those files are
already here.

## 2. Record your commit identity

From this directory:

```bash
git config user.name "YOUR NAME"
git config user.email "YOUR_GITHUB_EMAIL"
```

Use an email associated with your GitHub account, or your GitHub no-reply email.

## 3. Review and commit

```bash
git status --short
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[test]"
pytest -q
git add .
git diff --cached --check
git commit -m "Release flightcollapse v0.1.18 for collaborator evaluation"
git branch -M main
```

## 4. Connect and push

Replace `<OWNER>` with your GitHub username or organization:

```bash
git remote add origin https://github.com/<OWNER>/flightcollapse.git
git push -u origin main
```

GitHub may ask you to authenticate in a browser or with a personal access token.
Do not place a token in this repository or in a command saved to shell history.

## 5. Tag the evaluated version

```bash
git tag -a v0.1.18 -m "flightcollapse v0.1.18"
git push origin v0.1.18
```

Create a GitHub Release from the `v0.1.18` tag if you want a stable download
page. Use a title such as `flightcollapse v0.1.18 - collaborator evaluation` and
state that this is research software under evaluation.

## 6. Add collaborators

For a private repository, open **Settings -> Collaborators and teams**, invite
the evaluators, and give write access only to people who need to push branches.
Ask collaborators to install from the tag, not from a moving branch:

```bash
python -m pip install \
  "flightcollapse @ git+https://github.com/<OWNER>/flightcollapse.git@v0.1.18"
```

## Before making the repository public

- Confirm that the MIT license and author attribution are correct.
- Confirm that the four linked Claude artifacts are intended to remain
  link-accessible.
- Check that no patient data, read identifiers, BAMs, reference genomes,
  credentials, or private paths were added.
- Consider archiving the benchmark input tables or scripts separately if
  collaborators must reproduce the numerical claims.
