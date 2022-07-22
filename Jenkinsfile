pipeline {
  agent {
    docker {
      image 'node:lts'
      reuseNode true
    }

  }
  stages {
    stage('检出') {
      steps {
        checkout([
          $class: 'GitSCM',
          branches: [[name: GIT_BUILD_REF]],
          userRemoteConfigs: [[
            url: GIT_REPO_URL,
            credentialsId: CREDENTIALS_ID
          ]]])
        }
      }
      stage('同步到GitHub') {
        steps {
          sh 'git remote add github https://$GITHUB_NAME:$GITHUB_PASSWD@github.com/NekoMio/Blog'
          sh 'git push github HEAD:master'
        }
      }
      stage('拉取成品库') {
        steps {
          echo '从Coding拉取成品'
          script {
            sh "git clone https://${PROJECT_TOKEN_GK}:${PROJECT_TOKEN}@e.coding.net/WildRage/Blog/pages.git"
            sh "rm pages/* -rf"
          }

        }
      }
      stage('构建') {
        steps {
          echo '开始构建'
          sh 'pwd && ls'
          sh 'curl -fsSL https://bun.sh/install | bash'
          sh 'bun install'
          sh 'cd ./themes/suka/ && bun install --production'
          sh 'bun run build'
          sh 'cp ./_headers ./public/'
          sh 'cp ./public/* ./pages/ -rf'
        }
      }
      stage('Push To Coding') {
        steps {
          sh 'cd pages && git add .'
          sh 'cd pages && git commit -m "`date`"'
          sh 'cd pages && git push origin master'
        }
      }
      stage('Push To GitHub') {
        steps {
          sh 'cd pages && git remote add github https://$GITHUB_NAME:$GITHUB_PASSWD@github.com/NekoMio/nekomio.github.io'
          sh 'cd pages && git push github master'
        }
      }
    }
  }