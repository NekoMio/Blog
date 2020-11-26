pipeline {
  agent {
    docker {
      reuseNode true
      registryUrl 'https://coding-public-docker.pkg.coding.net'
      image 'public/docker/nodejs:14'
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
      stage('拉取成品库') {
        steps {
          echo '从Coding拉取成品'
          script {
            sh "git clone https://${PROJECT_TOKEN_GK}:${PROJECT_TOKEN}@e.coding.net/WildRage/Blog/pages.git"
            sh "cd pages"
            sh "rm * -rf"
            sh "cd .."
          }

        }
      }
      stage('构建') {
        steps {
          echo '开始构建'
          sh 'pwd && ls'
          sh 'npm install hexo-cli -g'
          sh 'yarn'
          sh 'hexo g'
          sh 'cp ./_headers ./public/ && cp ./aria2.html ./public/'
          sh 'cp ../public/* ./pages/ -rf'
        }
      }
      stage('Push To Coding') {
        steps {
          sh 'cd pages && git add .'
          sh 'git commit -m "`date`"'
          sh 'git push origin master'
        }
      }
    }
  }